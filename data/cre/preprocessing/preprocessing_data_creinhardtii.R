## Author: Francisco J. Romero-Campero, Ana Belén Romero-Losada
## Date: November 2018
## Contact: fran@us.es

## Documentation on how to generate org.db packages can be found on these links:
## http://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html
## 

cre.info <- read.table(file="Creinhardtii_281_v5.6.annotation_info.txt",header=T,comment.char = "",sep = "\t",as.is=T)
head(cre.info)

write.table(x = unique(cre.info$locusName),file = "cre_universe.txt",quote = F,row.names = F,col.names = F)


## Generate and write output data frame
id.cre.name <- data.frame(GID=cre.info$locusName,GENENAME=cre.info$locusName,stringsAsFactors = FALSE)
head(id.cre.name)

## Generate SYMBOL data frame
symbol.data.frame <- data.frame(GID=cre.info$locusName,SYMBOL=cre.info$locusName,stringsAsFactors = FALSE)
head(symbol.data.frame)
symbol.data.frame[1:20,]

## Remove duplicated rows
symbol.data.frame <- symbol.data.frame[!duplicated(symbol.data.frame),]
nrow(symbol.data.frame)

## Generate data frame for correspondence between Cre gene ids and Chlredraft gene ids used in KEGG
## Stablish correspondence betwee cre annotaiton and chlredraft
library(seqinr)
target.data <- read.fasta(file="GCF_000002595.1_v3.0_translated_cds.faa", seqtype="AA")

seq.annot <- getAnnot(target.data)
seq.names <- getName(target.data)

chlredraft <- vector(mode="character",length=length(seq.names))
for(i in 1:length(seq.annot))
{
  chlredraft[i] <- strsplit(strsplit(seq.annot[[i]],split="locXs_tag=")[[1]][2],split="]")[[1]][1]
}
names(chlredraft) <- seq.names

identity <- read.table(file="results_identity.txt",header = F,as.is = T)

get.first <- function(elto)
{
  return(elto[[1]])
}
cre.names <- sapply(strsplit(identity[[1]],split=".t"),get.first)

chlredraft.cre <- chlredraft[identity[[2]]]
names(chlredraft.cre) <- cre.names
sum(is.na(chlredraft.cre))

sorted.chlredraft <- chlredraft.cre[symbol.data.frame$GID]
sum(is.na(sorted.chlredraft))
names(sorted.chlredraft) <- NULL

## Generate CHLREDRAFT data frame
chlredraft.data.frame <- data.frame(GID=symbol.data.frame$GID,CHLREDRAFT=sorted.chlredraft,stringsAsFactors = FALSE)
head(chlredraft.data.frame)
chlredraft.data.frame[1:20,]

## Generate GO data.frame to create org.Db package
## Before reading the GO annotation file I have to replace the symbol ' by the word prime
## to avoid the error "EOF within quoted string"

## According to the vignettes of AnnotationForge bioconductor package for "Making Organisms
## Packages: "However to use the goTable argument, you have to follow a strict convention with 
## the data. Such a data.frame must have three columns only and these must correspond to the 
## gene id, GO id and evidence codes. These columns also have to be named as “GID”, “GO” and “EVIDENCE”

gid <- cre.info$locusName
go <- cre.info$GO
evidence <- rep("ISS",length(gid))

# for(i in 1:nrow(go.info))
# {
#  gid[i] <- gene.names[as.character(go.info$X.proteinId[i])]
#  go[i] <- go.info$goAcc[i]
# }


go.data.frame <- data.frame(GID=gid,GO=go,EVIDENCE=evidence,stringsAsFactors = FALSE)
head(go.data.frame)

go.data.frame <- subset(go.data.frame, GO != "")
head(go.data.frame)

go.expanded.data.frame <- c()
for(i in 1:nrow(go.data.frame))
{
  go.splitted <- strsplit(go.data.frame$GO[i],split=",")[[1]]

  if(length(go.splitted) > 0)
  {
    for(j in 1:length(go.splitted))
    {
      go.expanded.data.frame <- rbind(go.expanded.data.frame,c(go.data.frame[i,1],go.splitted[j]))
    }
  } else
  {
    go.expanded.data.frame <- rbind(go.expanded.data.frame,c(go.data.frame[i,1],""))
  }
}

head(go.expanded.data.frame)

go.data.frame <- data.frame(GID=go.expanded.data.frame[,1],
                            GO=go.expanded.data.frame[,2],
                            EVIDENCE=rep("ISS",nrow(go.expanded.data.frame)),stringsAsFactors = FALSE)
head(go.data.frame)


## Remove duplicated rows
go.data.frame <- go.data.frame[!duplicated(go.data.frame),]
go.data.frame[1:40,]
nrow(go.data.frame)
length(unique(go.data.frame$GID))

## For the rest of the data frames. According to the help info "the 1st column of every 
## data.frame must be labeled GID, and correspond to a gene ID that is universal for the 
## entire set of data.frames.  The GID is how the different tables will be joined internally"

## Generate enzyme data.frame

enzyme <- cre.info$ec
 
# for(i in 1:nrow(kegg.info))
# {
#   gid[i] <- gene.names[as.character(kegg.info$X.proteinId[i])]
#   enzyme[i] <- kegg.info$ecNum[i]
# }
 
enzyme.data.frame <- data.frame(GID=gid,ENZYME=enzyme,stringsAsFactors = FALSE)
head(enzyme.data.frame)

## Remove duplicated rows
enzyme.data.frame <- enzyme.data.frame[!duplicated(enzyme.data.frame),]
head(enzyme.data.frame)
nrow(enzyme.data.frame)
duplicated(enzyme.data.frame)
sum(duplicated(enzyme.data.frame))
length(unique(enzyme.data.frame$GID))


## Generate KOG data.frame
kog <- cre.info$KOG

# for(i in 1:nrow(kog.info))
# {
#   gid[i] <- gene.names[as.character(kog.info$proteinId[i])]
#   kog[i] <- kog.info$kogid[i]
# }

kog.data.frame <- data.frame(GID=gid,KOG=kog,stringsAsFactors = FALSE)
head(kog.data.frame)

## Remove duplicated rows
kog.data.frame <- kog.data.frame[!duplicated(kog.data.frame),]


## Generate KO data.frame
ko <- cre.info$KO

ko.data.frame <- data.frame(GID=gid,KO=ko,stringsAsFactors = FALSE)
head(ko.data.frame)

## Remove duplicated rows
ko.data.frame <- ko.data.frame[!duplicated(ko.data.frame),]
head(ko.data.frame)


## Generate PFAM data.frame
pfam <- cre.info$Pfam

pfam.data.frame <- data.frame(GID=gid,PFAM=pfam,stringsAsFactors = FALSE)
head(pfam.data.frame)

## Remove duplicated rows
pfam.data.frame <- pfam.data.frame[!duplicated(pfam.data.frame),]
head(pfam.data.frame)

## Generate PANTHER data.frame
panther <- cre.info$Panther

panther.data.frame <- data.frame(GID=gid,PANTHER=panther,stringsAsFactors = FALSE)
head(panther.data.frame)

## Remove duplicated rows
panther.data.frame <- panther.data.frame[!duplicated(panther.data.frame),]
head(panther.data.frame)

## Chlamydomonas reinhardtii Taxonomy ID: 3055

## Load require package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationForge", version = "3.8")
library(AnnotationForge)

makeOrgPackage(go=go.data.frame,
               SYMBOL=symbol.data.frame,
               ENZYME=enzyme.data.frame,
               KOG=kog.data.frame,
               KO=ko.data.frame,
               PFAM=pfam.data.frame,
               PANTHER=panther.data.frame,
               CHLREDRAFT=chlredraft.data.frame,
               version = "0.1",
               maintainer = "Francisco J. Romero-Campero <fran@us.es>",
               author = "Francisco J. Romero-Campero",
               outputDir = ".", 
               tax_id = "3055",
               genus = "Chlamydomonas",
               species = "reinhardtii",
               goTable = "go",
               verbose = TRUE)

install.packages("./org.Creinhardtii.eg.db/", repos=NULL)
remove.packages("org.Creinhardtii.eg.db")


library(org.Creinhardtii.eg.db)
columns(org.Creinhardtii.eg.db)
head(select(org.Creinhardtii.eg.db,columns = c("CHLREDRAFT"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID")))
head(select(org.Creinhardtii.eg.db,columns = c("ENZYME"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID")))
head(select(org.Creinhardtii.eg.db,columns = c("EVIDENCE"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID")))
head(select(org.Creinhardtii.eg.db,columns = c("EVIDENCEALL"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID")))
head(select(org.Creinhardtii.eg.db,columns = c("GO"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID")))
head(select(org.Creinhardtii.eg.db,columns = c("GOALL"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID")))
head(select(org.Creinhardtii.eg.db,columns = c("KOG"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID")))
head(select(org.Creinhardtii.eg.db,columns = c("ONTOLOGY"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID")))
head(select(org.Creinhardtii.eg.db,columns = c("ONTOLOGYALL"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID")))
head(select(org.Creinhardtii.eg.db,columns = c("SYMBOL"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID")))


library(org.At.tair.db)
columns(org.At.tair.db)
head(select(org.At.tair.db,columns = c("GO"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("ENZYME"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("PATH"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("GENENAME"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("SYMBOL"),keys=keys(org.At.tair.db,keytype = "TAIR")))


BiocManager::install("clusterProfiler", version = "3.8")
library(clusterProfiler)

cre.universe <- unique(select(org.Creinhardtii.eg.db,columns = c("GO"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID"))[["GID"]])
length(cre.universe)

cre.example <- read.table(file = "clusters/cre_trough_dark_light_peak_light_dark.txt",header = FALSE,as.is = TRUE)[[1]]

ego <- enrichGO(gene          = unique(c(cre.example)),
                universe      = cre.universe,
                OrgDb         = org.Creinhardtii.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE,
                keyType = "GID")

barplot(ego,drop=TRUE,showCategory = 10)

help(goplot)
goplot(x=ego,showCategory = 10)
dotplot(ego)
emapplot(ego)
cnetplot(ego)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview", version = "3.8")
library("pathview")

help(enrichKEGG)

cre.chlredraft.map <- select(org.Creinhardtii.eg.db,columns = c("CHLREDRAFT"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID"))
head(cre.chlredraft.map)
cre.ids <- cre.chlredraft.map$GID
chlredraft.ids <- cre.chlredraft.map$CHLREDRAFT
names(chlredraft.ids) <- cre.ids

cre.example.kegg <- chlredraft.ids[cre.example]
names(cre.example.kegg) <- NULL


kk <- enrichKEGG(gene = cre.example.kegg , organism = "cre",keyType = "kegg",
                 universe = cre.chlredraft.map$CHLREDRAFT,qvalueCutoff = 0.05)

draft.ids <- chlredraft.cre[cre.example]
names(draft.ids) <- NULL
draft.ids <- draft.ids[!is.na(draft.ids)]

draft.universe <- chlredraft.cre[cre.universe]
names(draft.universe) <- NULL
draft.universe <- draft.universe[!is.na(draft.universe)]

kk <- enrichKEGG(gene = draft.ids, organism = "cre",keyType = "kegg",
                 universe = draft.universe,qvalueCutoff = 0.05)





mkk <- enrichMKEGG(gene = draft.ids, organism = "cre",keyType = "kegg",universe = draft.universe,qvalueCutoff = 0.05)
head(mkk)

genes.pathway <- rep(0,length(draft.universe))
names(genes.pathway) <- draft.universe

genes.pathway[draft.ids] <- 1

res.kk <- head(as.data.frame(kk))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE), pathway.id =res.kk$ID[4], species = "cre",limit = list(gene=max(abs(genes.pathway)), cpd=1),gene.idtype ="kegg")


## Descargo de KEGG las proteinas anotadas con CHLREDRAFT 
## https://www.genome.jp/kegg-bin/show_organism?menu_type=genome_info&org=cre

## Preprocess gff3 to generate gtf with gene_id and transcript_id 
cre.gff3 <- read.table(file="../annotation/chlamydomonas_reinhardtii.gff3",header=F,quote = "#",as.is=T)
cre.gtf <- cre.gff3
head(cre.gff3)

unique(cre.gff3$V3)

for(i in 1:nrow(cre.gff3))
{
  current.attributes <- strsplit(cre.gff3$V9[i],split=";")[[1]]
  if(cre.gff3$V3[i] == "gene")
  {
    gene.id <- strsplit(current.attributes[2],split="=")[[1]][2]
    cre.gtf$V9[i] <- paste("gene_id", paste("\"",gene.id,"\";",sep=""))
  } else if(cre.gff3$V3[i] == "mRNA")
  {
    gene.id <- substr(strsplit(current.attributes[5],split="=")[[1]][2],start = 1,stop = 13)
    transcript.id <- strsplit(current.attributes[2],"=")[[1]][2]
    cre.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep="")), collapse = " ")
  } else if(cre.gff3$V3[i] == "exon")
  {
    gene.id <- substr(strsplit(current.attributes[1],split="=")[[1]][2],start = 1,stop = 13)
    trancript.id <- substr(strsplit(current.attributes[1],split="=")[[1]][2],start = 1,stop = 18)
    exon.number <- strsplit(strsplit(current.attributes[1],split="=")[[1]][2],split="exon.")[[1]][2]
    cre.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep=""),"exon_number",paste("\"",exon.number,"\";",sep="")), collapse = " ")
  } else if(cre.gff3$V3[i] == "five_prime_UTR" || cre.gff3$V3[i] == "three_prime_UTR" || cre.gff3$V3[i] == "CDS")
  {
    gene.id <- substr(strsplit(current.attributes[2],split="=")[[1]][2],start = 1,stop = 13)
    transcript.id <- substr(strsplit(current.attributes[2],split="=")[[1]][2],start = 1,stop = 18)
    cre.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep="")), collapse = " ")
  }
}

head(cre.gtf)

write.table(x = cre.gtf,file = "chlamydomonas_reinhardtii.gtf",sep = "\t",row.names = F,col.names = F,quote = F)

## Processing of gtf to remove common names
cre.gtf <- read.table(file="../annotation/chlamydomonas_reinhardtii.0.gtf",header=F,quote = "#",as.is=T,sep="\t")
cre.gtf.output <- cre.gtf

i <- 32

for(i in 1:nrow(cre.gtf))
{
  print(i)
  current.attributes <- strsplit(cre.gtf$V9[i],split=";")[[1]]
  if(cre.gtf$V3[i] == "mRNA")
  {
    gene.id <- substr(x=strsplit(current.attributes[2]," ")[[1]][3],start=2,stop=14)
    transcript.id <- strsplit(current.attributes[2]," ")[[1]][3]
    cre.gtf.output$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",
                                    transcript.id), collapse = " ")
  } 
}

write.table(x = cre.gtf.output,file = "chlamydomonas_reinhardtii.gtf",sep = "\t",row.names = F,col.names = F,quote = F)

## Generate TxDb package from gff3 file
library("GenomicFeatures")

## Generate chromosome.info data frame
library(seqinr)

cre.genome.data <- read.fasta(file = "../genome/chlamydomonas_reinhardtii.fa",seqtype = "DNA")
chromosome.names <- getName(cre.genome.data)
chromosome.lengths <- sapply(X=getSequence(cre.genome.data),FUN = length)

chromosomes.info <- data.frame(chrom=chromosome.names,length=chromosome.lengths,is_circular=FALSE)

## Meta data info
meta.data.info <- data.frame(name=c("Resource URL","Genome"),value=c("https://phytozome.jgi.doe.gov/","v5.0"))

cre.txdb <- makeTxDbFromGFF(file = "../annotation/chlamydomonas_reinhardtii.gff3",format = "gff3",dataSource = "Phytozome",organism = "Chlamydomonas reinhardtii",taxonomyId = 3055,chrominfo = chromosomes.info,metadata = meta.data.info)

cre.txdb
genes(cre.txdb)

?makeTxDbPackage

makeTxDbPackage(txdb = cre.txdb, version = "0.1", maintainer = "Francisco J. Romero-Campero <fran@us.es>", author = "Francisco J. Romero-Campero")

install.packages("./TxDb.Creinhardtii.Phytozome/", repos=NULL)
## loading packages
library(ChIPseeker)
library(TxDb.Creinhardtii.Phytozome)
txdb <- TxDb.Creinhardtii.Phytozome
library(clusterProfiler)

files <- c("peaks_H3K27me3_before_starvation_1_peaks.narrowPeak","peaks_H3K27me3_sulfur_starvation_1_peaks.narrowPeak")
peak <- readPeakFile(files[[2]])
peak

covplot(peak, weightCol="X109")

covplot(peak, weightCol="X109",chrs="chromosome_1")
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
tagMatrix <- getTagMatrix(peak, windows=promoter)

tagHeatmap(tagMatrix, xlim=c(-1000, 1000), color="red")

plotAvgProf(tagMatrix, xlim=c(-1000, 1000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf(tagMatrix, xlim=c(-1000, 1000), conf = 0.95, resample = 1000)

peakAnno <- annotatePeak(files[[2]], tssRegion=c(-1000, 1000),
                         TxDb=txdb, annoDb="org.Creinhardtii.eg.db")

plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

peak.annotation <- as.data.frame(peakAnno)
head(peak.annotation)
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-1000, 1000))

plotAvgProf(tagMatrixList, xlim=c(-1000, 1000), conf=0.95,resample=500, facet="row")

tagHeatmap(tagMatrixList, xlim=c(-1000, 1000), color=NULL)