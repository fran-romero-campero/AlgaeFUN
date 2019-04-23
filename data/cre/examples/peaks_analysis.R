## Author: Ana Belén Romero-Losada,  Francisco J. Romero-Campero
## Date: April 2019
## Contact: fran@us.es

## Loading packages
library(ChIPseeker)
library(ChIPpeakAnno)
library(rtracklayer)

library(TxDb.Creinhardtii.Phytozome)

library(clusterProfiler)

## Auxiiary functions
extract.annotation <- function(annotation.str)
{
  return(strsplit(x = annotation.str,split=" \\(")[[1]][1])
}

input <- list(analysis = "genomic_locations", 
              peaks="output_peaks.narrowPeak", 
              promoter_length=1000,
              tes_length=1000,
              selected_genomic_features=c("Promoter","Gene Body"),
              bw_file="chip.bw")
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
              title="Distribution of genomic loci relative to TSS",ylab = "Genomic Loci (%) (5' -> 3')")
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

## Metagene plot

## Extraction of the genomic features of the specified genes.
genes.data <- subset(genes(txdb), gene_id %in% genes)

## Extraction of the TSS 
genes.tss <- resize(genes.data, width=1, fix='start')

## Centering around TSS with promoter length
around.genes.tss <- genes.tss
start(around.genes.tss) <- start(genes.tss) - input$promoter_length
end(around.genes.tss) <- end(genes.tss) + input$promoter_length
feature.recentered <- around.genes.tss

## Importing bigWig file
cvglists <- sapply(input$bw_file, import, 
                   format="BigWig", 
                   which=feature.recentered, 
                   as="RleList")

## Extracting the signal around TSS with promoter length
number.tiles <- 2*input$promoter_length/10
tss.sig <- featureAlignedSignal(cvglists, around.genes.tss, 
                                upstream=input$promoter_length, 
                                downstream=input$promoter_length,
                                n.tile=number.tiles) 

## Plotting the results around TSS (you may have to change the names of the conditions)
profile.to.plot <- colMeans(tss.sig[[1]],na.rm = TRUE)
plot(profile.to.plot,type="l",col="blue",lwd=3,ylab="",cex.lab=2,axes=FALSE,xlab="")#,ylim=c(0,830))
polygon(c(1,1:length(profile.to.plot),length(profile.to.plot)),
        c(0,profile.to.plot,0),col="lightblue")

axis(side = 1,
     labels = c(-input$promoter_length,-input$promoter_length/2,
                "TSS",
                input$promoter_length/2,input$promoter_length),
     at = c(1,number.tiles/4,number.tiles/2,3*number.tiles/4,number.tiles),lwd=2,cex=1.5,las=2,cex=2)

## Extraction of the TES
genes.tes <- resize(genes.data, width=1, fix='end')
head(genes.tes)

## Centering around TES 
around.genes.tes <- genes.tes
start(around.genes.tes) <- start(genes.tes) - input$tes_length
end(around.genes.tes) <- end(genes.tes) + input$tes_length

## Extracting the signal 2kb from the center with 50 tile (This may be slow)
number.tiles.tes <- 2 * input$tes_length /10
tes.sig <- featureAlignedSignal(cvglists, around.genes.tes, 
                                upstream=input$tes_length, 
                                downstream=input$tes_length,
                                n.tile=number.tiles.tes) 

## Plotting the results around TES
profile.to.plot <- colMeans(tes.sig[[1]],na.rm = TRUE)
plot(profile.to.plot,type="l",col="blue",lwd=3,ylab="",cex.lab=2,axes=FALSE,xlab="")#,ylim=c(0,830))
polygon(c(1,1:length(profile.to.plot),length(profile.to.plot)),
        c(0,profile.to.plot,0),col="lightblue")

axis(side = 1,
     labels = c(-input$tes_length,-input$tes_length/2,
                "TES",
                input$tes_length/2,input$tes_length),
     at = c(1,number.tiles.tes/4,number.tiles.tes/2,3*number.tiles.tes/4,number.tiles.tes),lwd=2,cex=1.5,las=2,cex=2)


## Extracting signal for the entire gene
genes <- genes[1:1000]

gene.body.signal <- matrix(0,nrow=length(genes),ncol=100)

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
    sig <- featureAlignedSignal(cvglists, current.gene.data, #n.tile=100)
                                upstream=recenter.value[1], downstream=recenter.value[1],n.tile=100) 
    
    
    gene.body.signal[i,] <- sig[[1]][1,]
  }
}

## Computing the average signal for each condition
mean.gene.body.signal <- colMeans(gene.body.signal)

## Merging the signal from the TSS, the gene body and TES
tss.profile <- colMeans(tss.sig[[1]],na.rm = TRUE)
tes.profile <- colMeans(tes.sig[[1]],na.rm = TRUE)


tss.profile.to.plot <- tss.profile[1:(length(tss.profile)/2)]
tes.profile.to.plot <- tes.profile[1:(length(tes.profile)/2)]

length(tss.profile.to.plot)
length(tes.profile.to.plot)

n <- 4


step.to.take.tss <- length(tss.profile.to.plot)/(length(mean.gene.body.signal)/n)
step.to.take.tes <- length(tes.profile.to.plot)/(length(mean.gene.body.signal)/n)

tss.profile.to.plot <- tss.profile.to.plot[seq(from=1,to=length(tss.profile),by=step.to.take.tss)]
tes.profile.to.plot <- tes.profile.to.plot[seq(from=1,to=length(tes.profile),by=step.to.take.tes)]

mean.merged <- c(tss.profile.to.plot,
                 mean.gene.body.signal,
                 tes.profile.to.plot)

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



mut.merged <- c(colMeans(signal.around.tss[["clfswn"]],na.rm = TRUE)[1:25],mean.mut,colMeans(signal.around.tes[["clfswn"]],na.rm = TRUE)[26:50])



