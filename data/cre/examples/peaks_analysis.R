## Author: Ana Belén Romero-Losada,  Francisco J. Romero-Campero
## Date: April 2019
## Contact: fran@us.es

## Loading packages
library(ChIPseeker)
library(ChIPpeakAnno)
library(rtracklayer)
library(seqinr)

library(TxDb.Creinhardtii.Phytozome)

library(clusterProfiler)

## Auxiiary functions
extract.annotation <- function(annotation.str)
{
  return(strsplit(x = annotation.str,split=" \\(")[[1]][1])
}

## Function to compute the reverse complement
reverse.complement <- function(dna.sequence)
{
  return(c2s(comp(rev(s2c(dna.sequence)),forceToLower = FALSE)))
}


input <- list(analysis = "genomic_locations", 
              peaks="output_peaks.narrowPeak", 
              promoter_length=1000,
              tes_length=1000,
              selected_genomic_features=c("Promoter","Gene Body"),
              bw_file="chip.bw",
              min_score_pwm = 95)
# gene_set

## Define transcript data base for the corresponding microalgae
txdb <- TxDb.Creinhardtii.Phytozome

## Load peak file
peaks <- readPeakFile(peakfile = input$peaks,header=FALSE)

## Define promoter region around TSS
promoter <- getPromoters(TxDb=txdb, upstream=input$promoter_length, downstream=input$promoter_length)

## Plot signal in promoter
#plotAvgProf(tagMatrix, xlim=c(-input$promoter_length, input$promoter_length),
#            xlab="Genomic Region (5'->3')", ylab = "Signal")



## Annotate genomic loci
peakAnno <- annotatePeak(peaks, tssRegion=c(-input$promoter_length, input$promoter_length),
                         TxDb=txdb)#, annoDb="org.Creinhardtii.eg.db")

## Plot pie chart with annotation
plotAnnoPie(peakAnno)
#plotAnnoBar(peakAnno)
#vennpie(peakAnno)
#upsetplot(peakAnno)
#upsetplot(peakAnno, vennpie=TRUE)

## Plot distance to tss
plotDistToTSS(peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")
dev.off()
## Extract genes 
peak.annotation <- as.data.frame(peakAnno)
simple.annotation <- sapply(X = as.vector(peak.annotation$annotation),FUN = extract.annotation)
names(simple.annotation) <- NULL

genes.promoter <- peak.annotation$geneId[simple.annotation == "Promoter"]
genes.5utr <- peak.annotation$geneId[simple.annotation == "5' UTR"]
genes.3utr <- peak.annotation$geneId[simple.annotation == "3' UTR"]
genes.exon <- peak.annotation$geneId[simple.annotation == "Exon"]
genes.intron <- peak.annotation$geneId[simple.annotation == "Intron"]

## Select final gene set
genes <- c()
if( "Promoter" %in% input$selected_genomic_features )
{
  genes <- c(genes,genes.promoter)
}

if("Gene Body" %in% input$selected_genomic_features)
{
  genes <- c(genes,genes.5utr,genes.3utr,genes.exon,genes.intron)
}

genes <- unique(genes)

length(genes)

## Extraction of the genomic features of the specified genes.
genes.data <- subset(genes(txdb), gene_id %in% genes)
genes.data


## Metagene plot

## Extraction of the TSS 
genes.tss <- resize(genes.data, width=1, fix='start')
head(genes.tss)
## Centering around TSS with promoter length
around.genes.tss <- genes.tss
start(around.genes.tss) <- start(genes.tss) - input$promoter_length
end(around.genes.tss) <- end(genes.tss) + input$promoter_length

## Importing bigWig file
cvglists <- sapply(input$bw_file, import, 
                   format="BigWig", 
                   which=around.genes.tss, 
                   as="RleList")

## Extracting the signal around TSS with promoter length
number.tiles <- 2*input$promoter_length/20
tss.sig <- featureAlignedSignal(cvglists, around.genes.tss, 
                                upstream=input$promoter_length, 
                                downstream=input$promoter_length,
                                n.tile=number.tiles) 

## Plotting the results around TSS (you may have to change the names of the conditions)
profile.around.tss <- colMeans(tss.sig[[1]],na.rm = TRUE)
plot(profile.around.tss,type="l",col="blue",lwd=3,ylab="",cex.lab=2,axes=FALSE,xlab="")
polygon(c(1,1:length(profile.around.tss),length(profile.around.tss)),
        c(0,profile.around.tss,0),col="lightblue")

axis(side = 1,
     labels = c(-input$promoter_length,-input$promoter_length/2,
                "TSS",
                input$promoter_length/2,input$promoter_length),
     at = c(1,number.tiles/4,number.tiles/2,3*number.tiles/4,number.tiles),lwd=2,cex=1.5,las=2,cex=2)

## Extraction of the TES
genes.tes <- resize(genes.data, width=1, fix='end')

## Centering around TES 
around.genes.tes <- genes.tes
start(around.genes.tes) <- start(genes.tes) - input$tes_length
end(around.genes.tes) <- end(genes.tes) + input$tes_length

## Extracting the signal 2kb from the center with 50 tile (This may be slow)
number.tiles.tes <- 2 * input$tes_length /20
tes.sig <- featureAlignedSignal(cvglists, around.genes.tes, 
                                upstream=input$tes_length, 
                                downstream=input$tes_length,
                                n.tile=number.tiles.tes) 

## Plotting the results around TES
profile.around.tes <- colMeans(tes.sig[[1]],na.rm = TRUE)
plot(profile.around.tes,type="l",col="red4",lwd=3,ylab="",cex.lab=2,axes=FALSE,xlab="")#,ylim=c(0,830))
polygon(c(1,1:length(profile.around.tes),length(profile.around.tes)),
        c(0,profile.around.tes,0),col="bisque")

axis(side = 1,
     labels = c(-input$tes_length,-input$tes_length/2,
                "TES",
                input$tes_length/2,input$tes_length),
     at = c(1,number.tiles.tes/4,number.tiles.tes/2,3*number.tiles.tes/4,number.tiles.tes),lwd=2,cex=1.5,las=2,cex=2)


## Extracting signal for the entire gene
genes <- genes[1:100]

gene.body.signal <- matrix(0,nrow=length(genes),ncol=100)
i <- 1
for(i in 1:length(genes))
{
  ## Printing current gene and extracting the features for the next two 
  ## genes
  print(i)
  current.gene.data <- genes.data[c(i,i+1),]
  
  ## Values for recentering (length of each feature divided by two)
  recenter.value <- floor(width(current.gene.data)/2)
  
  ## Centering features
  centered.current.gene.data <- current.gene.data
  start(centered.current.gene.data[1,]) <- start(centered.current.gene.data[1,]) + recenter.value[1]
  start(centered.current.gene.data[2,]) <- start(centered.current.gene.data[2,]) + recenter.value[2]
  end(centered.current.gene.data) <- start(centered.current.gene.data)
  
  start(current.gene.data[2,]) <- start(centered.current.gene.data[2,]) - recenter.value[1]
  end(current.gene.data[2,]) <- end(centered.current.gene.data[2,]) + recenter.value[1] -1
  
  ## Getting widths of the features
  widths <- width(current.gene.data)
  
  ## Getting the signal for the current gene data
  cvglists <- sapply(input$bw_file, import, 
                     format="BigWig", 
                     which=current.gene.data, 
                     as="RleList")
  ## You may want to change the names of the conditions here
  
  #names(cvglists) <- c("Col-0", "atbmi1abc") # example atbmi1abc vs Col-0
  if(widths[1] > 100)
  {
    sig <- featureAlignedSignal(cvglists, 
                                current.gene.data, #n.tile=100)
                                upstream=recenter.value[1], 
                                downstream=recenter.value[1],
                                n.tile=100) 
    
    
    gene.body.signal[i,] <- sig[[1]][1,]
  }
}

## Computing the average signal for each condition
mean.gene.body.signal <- colMeans(gene.body.signal)

profile.around.tes

## Merging the signal from the TSS, the gene body and TES
tss.profile <- colMeans(tss.sig[[1]],na.rm = TRUE)
tes.profile <- colMeans(tes.sig[[1]],na.rm = TRUE)

length(tes.profile)


tes.profile[100]
tes.profile[101]
mean.gene.body.signal[100]


tss.profile.to.plot <- tss.profile[1:(length(tss.profile)/2)]
tes.profile.to.plot <- tes.profile[((length(tes.profile)/2)+1):length(tes.profile)]

length(tss.profile.to.plot)
length(tes.profile.to.plot)

n <- 4


step.to.take.tss <- length(tss.profile.to.plot)/(length(mean.gene.body.signal)/n)
step.to.take.tes <- length(tes.profile.to.plot)/(length(mean.gene.body.signal)/n)

tss.profile.to.plot.2 <- tss.profile.to.plot[seq(from=1,to=length(tss.profile.to.plot),by=step.to.take.tss)]
tes.profile.to.plot.2 <- tes.profile.to.plot[seq(from=1,to=length(tes.profile),by=step.to.take.tes)]

mean.merged <- c(tss.profile.to.plot.2,
                 mean.gene.body.signal,
                 tes.profile.to.plot.2)

plot(mean.merged,type="l",col="blue",lwd=3,ylab="",cex.lab=2,axes=FALSE,xlab="")#,ylim=c(0,830))
polygon(c(1,1:length(mean.merged),length(mean.merged)),
        c(0,mean.merged,0),col="lightblue")

axis(side = 1,
     labels = c(-input$tes_length,
                -input$tes_length/2,
                "TES",
                input$tes_length/2,
                input$tes_length),
     at = c(1,
            number.tiles.tes/4,
            number.tiles.tes/2,
            3*number.tiles.tes/4,
            number.tiles.tes),
     lwd=2,cex=1.5,las=2,cex=2)



## Visualization of individual genes
genes.data <- as.data.frame(genes.data)
exons.data <- as.data.frame(exons(txdb))
cds.data <- as.data.frame(cds(txdb))

gene.name <- genes[1]

target.gene.body <- genes.data[gene.name,]
target.gene.chr <- as.character(target.gene.body$seqnames)
target.gene.start <- target.gene.body$start
target.gene.end <- target.gene.body$end

target.gene.strand <- as.character(target.gene.body$strand)

## Extract cds annotation
cds.data.target.gene <- subset(cds.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))

## Extract exons annotation
exons.data.target.gene <- subset(exons.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))

## Determine the genome range to plot including promoter, gene body and 5' UTR
## This depends on whether the gene is on the forward or reverse strand
range.to.plot <- target.gene.body

if(target.gene.strand == "+")
{
  range.to.plot$start <- range.to.plot$start - input$promoter_length
  range.to.plot$end <- range.to.plot$end + input$tes_length
} else if (target.gene.strand == "-")
{
  range.to.plot$end <- range.to.plot$end + input$promoter_length
  range.to.plot$start <- range.to.plot$start - input$tes_length
}

## Compute the length of the genome range to represent
current.length <- range.to.plot$end - range.to.plot$start

## Compute profile in gene
selected.bigwig.files <- input$bw_file
selected.bed.files <- input$peaks

## Since ChIPpeakAnno needs more than one region to plot our region
## is duplicated 
regions.plot <- GRanges(rbind(range.to.plot,range.to.plot))

## Import signal from the bigwig files
cvglists <- sapply(selected.bigwig.files, import, 
                   format="BigWig", 
                   which=regions.plot, 
                   as="RleList")

## Compute signal in the region to plot
chip.signal <- featureAlignedSignal(cvglists, regions.plot, 
                                    upstream=ceiling(current.length/2), 
                                    downstream=ceiling(current.length/2),
                                    n.tile=current.length) 

## Compute mean signal 
if(target.gene.strand == "+")
{
  chip.signal.mean <- colMeans(chip.signal[[1]],na.rm = TRUE)
} else if (target.gene.strand == "-")
{
  chip.signal.mean <- rev(colMeans(chip.signal[[1]],na.rm = TRUE))
}

## Normalization
chip.signal.mean <- 20 * chip.signal.mean / max(chip.signal.mean)

## Determine upper limit of the graph
upper.lim <- 21

## Draw DNA strand
gene.height <- -25
cord.x <- 1:current.length

plot(cord.x, rep(gene.height,length(cord.x)),type="l",col="black",lwd=3,ylab="",
     cex.lab=2,axes=FALSE,xlab="",main="",cex.main=2,
     ylim=c(-35,upper.lim),xlim=c(-3000,max(cord.x)))
  
## Extract exons for target gene
exons.data.target.gene <- subset(exons.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
  
## Transform exon coordinates to current range
min.pos <- min(exons.data.target.gene$start)
  
if(target.gene.strand == "+")
{
  exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + input$promoter_length
  exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + input$promoter_length
} else if(target.gene.strand == "-")
{
  exons.data.target.gene$start <- exons.data.target.gene$start - min.pos + input$tes_length
  exons.data.target.gene$end <- exons.data.target.gene$end - min.pos + input$tes_length
}
  
## Represent exons
exon.width <- 2
for(i in 1:nrow(exons.data.target.gene))
{
  # Determine start/end for each exon
  current.exon.start <- exons.data.target.gene$start[i]
  current.exon.end <- exons.data.target.gene$end[i]
  
  ## Determine coordinates for each exon polygon and represent it
  exon.x <- c(current.exon.start,current.exon.end,current.exon.end,current.exon.start)
  exon.y <- c(gene.height + exon.width, gene.height + exon.width, gene.height - exon.width, gene.height - exon.width)
  
  polygon(x = exon.x, y = exon.y, col = "blue",border = "blue")
}
  
## Extract cds for target gene
cds.data.target.gene <- subset(cds.data, seqnames == target.gene.chr & (start >= target.gene.start & end <= target.gene.end))
  
## Transform cds coordinates to current range
if(target.gene.strand == "+")
{
  cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + input$promoter_length
  cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + input$promoter_length
} else if (target.gene.strand == "-")
{
  cds.data.target.gene$start <- cds.data.target.gene$start - min.pos + input$tes_length
  cds.data.target.gene$end <- cds.data.target.gene$end - min.pos + input$tes_length
}
  
cds.width <- 3
for(i in 1:nrow(cds.data.target.gene))
{
  # Determine current cds start/end
  current.cds.start <- cds.data.target.gene$start[i]
  current.cds.end <- cds.data.target.gene$end[i]
  
  # Determine curret cds coordinates for the polygon and represent it
  cds.x <- c(current.cds.start,current.cds.end,current.cds.end,current.cds.start)
  cds.y <- c(gene.height + cds.width, gene.height + cds.width, gene.height - cds.width, gene.height - cds.width)
  
  polygon(x = cds.x, y = cds.y, col = "blue",border = "blue")
}
  
## Draw arrow to represent transcription direction 
if(target.gene.strand == "+")
{
  lines(c(input$promoter_length,input$promoter_length,input$promoter_length+100),y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
  lines(c(input$promoter_length+50,input$promoter_length+100),y=c(gene.height+6,gene.height+5),lwd=3)
  lines(c(input$promoter_length+50,input$promoter_length+100),y=c(gene.height+4,gene.height+5),lwd=3)
} else if (target.gene.strand == "-")
{
  lines(c(current.length - input$promoter_length, current.length - input$promoter_length, current.length - input$promoter_length-100),y=c(gene.height,gene.height+5,gene.height+5),lwd=3)
  lines(c(current.length - input$promoter_length-50, current.length - input$promoter_length - 100),y=c(gene.height + 6, gene.height + 5),lwd=3)
  lines(c(current.length - input$promoter_length-50, current.length - input$promoter_length - 100),y=c(gene.height + 4, gene.height + 5),lwd=3)
}
  
## Draw promoter range
if(target.gene.strand == "+")
{
  axis(side = 1,labels = c(- input$promoter_length, - input$promoter_length / 2,"TSS"),at = c(1,input$promoter_length/2,input$promoter_length),lwd=2,cex=1.5,las=2,cex=2)
} else if(target.gene.strand == "-")
{
  axis(side = 1,labels = c("TSS",- input$promoter_length / 2,- input$promoter_length),at = c(current.length-input$promoter_length,current.length-input$promoter_length/2, current.length),lwd=2,cex=1.5,las=2,cex=2)
}
  
## Draw gene name

text(x = current.length / 2, y = -33 , 
     labels = bquote(italic(.(gene.name))),cex = 1.7,font = 3)

## Draw peak regions
## Width of the rectangule representing the peak reagion
peak.width <- 1

## Extract bed file name 1 and read it
current.peaks <- read.table(file=input$peaks,header = F, as.is = T)
peak.coordinates <- subset(current.peaks, V1 == range.to.plot$seqnames & V2 >= range.to.plot$start & V3 <= range.to.plot$end) 
current.peaks.to.plot <- peak.coordinates[,2:3]

## Transform coordinates 
current.peaks.to.plot <- current.peaks.to.plot - range.to.plot$start
  
## Check if there are peaks for the target gene
if(nrow(current.peaks.to.plot) > 0)
{
  #motifs.in.peaks <- vector(mode="list", length=nrow(current.peaks.to.plot))
  for(j in 1:nrow(current.peaks.to.plot))
  {
    ## Extract start and end point of each peak region
    current.peak.start <- current.peaks.to.plot[j,1]
    current.peak.end <- current.peaks.to.plot[j,2]
      
    ## Computer coordinates for polygon and draw it
    peak.x <- c(current.peak.start,current.peak.end,
                current.peak.end,current.peak.start)
    peak.y <- c(peak.width - 12,   peak.width - 12, 
              - peak.width - 12, - peak.width - 12)  

    polygon(x = peak.x, y = peak.y, col = "bisque", border = "darkred",lwd=2)
  }
}

line.colors <- "blue"
area.colors <- "lightblue"

## Draw profiles
## Compute base line for current TF
current.base.line <- - 10
  
## Represent signal from the current TF
lines(chip.signal.mean+current.base.line,type="l",col=line.colors,lwd=3)
  
## Determine polygon coordinates and represent it
cord.y <- c(current.base.line,chip.signal.mean+current.base.line,current.base.line)
cord.x <- 1:length(cord.y)
  
polygon(cord.x,cord.y,col=area.colors)

## Load microalgae genome
microalgae.genome.data <- read.fasta(file = "../genome/chlamydomonas_reinhardtii.fa",seqtype = "DNA")
microalgae.genome <- getSequence(microalgae.genome.data)
names(microalgae.genome) <- getName(microalgae.genome.data)

## Load Position Weight Matrices
## Open file connection
con <- file("../../../web_app/jaspar_motifs/pfm_plants_20180911.txt",open = "r")
  
## Empty list for storing PWM
motifs.pwm <- vector(mode="list",length = 453)
motif.ids <- vector(mode="character",length=453)
motif.names <- vector(mode="character",length=453)
  
## Load 64 PWM
for(j in 1:453)
{
  ## First line contains motif id and name
  first.line <- readLines(con,1)
    
  motif.ids[j] <- strsplit(first.line,split=" ")[[1]][1]
  motif.names[j] <- strsplit(first.line,split=" ")[[1]][2]
    
  ## Next four line contians probabilites for each nucleotide
  a.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  c.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  g.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
  t.row <- as.numeric(strsplit(readLines(con,1),split="( )+")[[1]])
    
  ## Construct PWM
  motif.pwm <- matrix(nrow = 4,ncol=length(a.row))
  
  motif.pwm[1,] <- a.row
  motif.pwm[2,] <- c.row 
  motif.pwm[3,] <- g.row
  motif.pwm[4,] <- t.row
    
  rownames(motif.pwm) <- c("A","C","G","T")
  
  motifs.pwm[[j]] <- prop.table(motif.pwm,2)
}
  
## Close file connection
close(con)
  
## Naming list with PWM
names(motifs.pwm) <- motif.names
names(motif.ids) <- motif.names
  
## Draw peak regions for each TF and determing TF binding sequences
  
## Determine TFBS motifs to search for
if(input$all.motifs)
{
  selected.motifs.pwm <- motifs.pwm
} else
{
  selected.motifs.pwm <- motifs.pwm[input$selected.motifs]
}
  
selected.motif.names <- names(selected.motifs.pwm)
selected.motif.ids <- motif.ids[selected.motif.names]
  
## Initialize data frame containing TF binding sequences in the peak regions
df.hits <- data.frame(0,"","","")
colnames(df.hits) <- c("position","id","name","seq")
  
## Identify TF binding DNA motifs 
if(nrow(current.peaks.to.plot) > 0)
{
  #motifs.in.peaks <- vector(mode="list", length=nrow(current.peaks.to.plot))
  for(j in 1:nrow(current.peaks.to.plot))
  {
    ## Genomic coordinates of the current peak
    peak.chr <- peak.coordinates[j, 1]
    peak.start <- peak.coordinates[j, 2]
    peak.end <- peak.coordinates[j, 3]
    
    ## Extract start and end point of each peak region in our plot
    current.peak.start <- current.peaks.to.plot[j,1]
    current.peak.end <- current.peaks.to.plot[j,2]
    
        
    ## Extract peak sequence
    peak.sequence <- c2s(microalgae.genome[[peak.chr]][peak.start:peak.end])
    peak.rev.comp.sequence <- reverse.complement(peak.sequence)
      
    for(k in 1:length(selected.motifs.pwm))
    {
      print(k)
      motif.pwm <- selected.motifs.pwm[[k]]
          
      hits.fw <- matchPWM(motif.pwm, peak.sequence, 
                          min.score = paste0(input$min_score_pwm,"%"))
      hits.fw
      hits.fw.seqs <- as.data.frame(hits.fw)[[1]]
      hits.fw <- as(hits.fw, "IRanges")
      hits.fw.start <- start(hits.fw)
      hits.fw.end <- end(hits.fw)
          
      if(length(hits.fw.start) > 0)
      {
        df.hits.fw <- data.frame(((hits.fw.start+hits.fw.end)/2) + current.peak.start,
                                 rep(selected.motif.ids[k],length(hits.fw.start)),
                                 rep(selected.motif.names[k],length(hits.fw.start)),
                                 hits.fw.seqs)
        colnames(df.hits.fw)  <- c("position","id","name","seq")
        df.hits <- rbind(df.hits,df.hits.fw)
      }
          
      hits.rev <- matchPWM(motif.pwm, peak.rev.comp.sequence, 
                           min.score = paste0(input$min.score.pwm,"%"))
      hits.rev.seqs <- as.data.frame(hits.rev)[[1]]
      hits.rev.seqs <- sapply(hits.rev.seqs,reverse.complement)
      names(hits.rev.seqs) <- NULL
          
      hits.rev <- as(hits.rev, "IRanges")
      hits.rev.start <- nchar(peak.sequence) - end(hits.rev) + 1
      hits.rev.end <- nchar(peak.sequence) - start(hits.rev) + 1
          
      if(length(hits.rev.start) > 0)
      {
        df.hits.rev <- data.frame(((hits.rev.start+hits.rev.end)/2) + current.peak.start,
                                  rep(selected.motif.ids[k],length(hits.rev.start)),
                                  rep(selected.motif.names[k],length(hits.rev.start)),
                                  hits.rev.seqs)
        colnames(df.hits.rev)  <- c("position","id","name","seq")
        df.hits <- rbind(df.hits,df.hits.rev)
      }
    }
  }
}

## Remove first line of the data frame added just for technical reason
df.hits <- df.hits[-1,]
nrow(df.hits)
  
## Draw TF binding sites
detected.tfbs <- unique(as.vector(df.hits$name))

## TF binding sites colors and symbol shapes
symbol.shapes <- c(17, 18, 19, 15)
symbol.color <- c("blue", "red", "darkgreen", "magenta")
  
number.of.shapes <- ceiling(length(detected.tfbs) / length(symbol.color))
necessary.shapes <- rep(symbol.shapes[1:number.of.shapes],each = length(detected.tfbs)/number.of.shapes)
necessary.colors <- rep(symbol.color,number.of.shapes)
  
if(length(detected.tfbs) > 0)
{
  for(i in 1:length(detected.tfbs))
  {
    current.tfbs <- detected.tfbs[i]
    current.shape <- necessary.shapes[i]
    current.color <- necessary.colors[i]
      
    positions <- subset(df.hits, name == current.tfbs)
      
    for(j in 1:nrow(positions))
    {
      pos.to.draw <- positions$position[j]
      
      points(x = pos.to.draw, y = -15,
             pch = current.shape, col = current.color, cex = 1)
    }
  }
    
  ## Add legend for TFBS
  legend.step <- 5
  for(i in 1:length(detected.tfbs))
  {
    points(x = -3000, y = upper.lim - (i-1)*legend.step, 
           pch=necessary.shapes[i], col = necessary.colors[i],cex = 1)
      
      
    current.seq <- as.character(subset(df.hits,name == detected.tfbs[i])[["seq"]][[1]])
    current.label <- paste(c(detected.tfbs[i], "  -  ", current.seq ),collapse="")
      
    text(x = -2900, y = upper.lim - (i-1)*legend.step, labels = current.label,
         adj = 0,cex = 0.7)
  }
}
