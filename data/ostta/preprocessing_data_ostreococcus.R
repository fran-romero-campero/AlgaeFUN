## Author: Francisco J. Romero-Campero, Ana Belén Romero-Losada
## Date: November 2018
## Contact: fran@us.es

## Documentation on how to generate org.db packages can be found on these links:
## http://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html
## 

## Establish correspondence between protein ids and ostta gene/protein names
library("seqinr")

## Read fasta file and split sequence names with the token |
proteins.info <- read.fasta(file = "Ostta4221_3_GeneCatalog_proteins_20161028.aa.fasta",seqtype = "AA")
seq.names <- getName(proteins.info)
splitted.seq.names <- strsplit(seq.names,split="\\|")

## Loop to extract protein ids and ostta names. Initialize accumulators
proteins.id <- vector(mode="character",length=length(seq.names))
ostta.names <- vector(mode="character",length=length(seq.names))

for(i in 1:length(seq.names))
{
  proteins.id[i] <- splitted.seq.names[[i]][3]
  ostta.names[i] <- substr(x = splitted.seq.names[[i]][4],start = 4,stop = 16)
}

## Generate and write output data frame
id.ostta.name <- data.frame(GID=proteins.id,GENENAME=ostta.names,SYMBOL=ostta.names)
head(id.ostta.name)

write.table(x = id.ostta.name,file = "correspondence_gene_id_ostta_name.tsv",sep="\t",quote=F,row.names = F)

id.ostta.name <- read.table(file = "correspondence_gene_id_ostta_name.tsv",sep="\t",header=T,as.is=T)
head(id.ostta.name)

gene.ids <- id.ostta.name$id
gene.names <- id.ostta.name$name
names(gene.names) <- gene.ids

## Generate GO data.frame to create org.Db package
## Before reading the GO annotation file I have to replace the symbol ' by the word prime
## to avoid the error "EOF within quoted string"
go.info <- read.table(file="Ostta4221_3_GeneCatalog_proteins_20161028_GO.tab",header = T,as.is = T,sep = "\t",comment.char = "")
head(go.info)

## According to the vignettes of AnnotationForge bioconductor package for "Making Organisms
## Packages: "However to use the goTable argument, you have to follow a strict convention with 
## the data. Such a data.frame must have three columns only and these must correspond to the 
## gene id, GO id and evidence codes. These columns also have to be named as “GID”, “GO” and “EVIDENCE”

gid <- vector(mode="character",length=nrow(go.info))
go <- vector(mode="character",length=nrow(go.info))
evidence <- rep("ISS",nrow(go.info))

for(i in 1:nrow(go.info))
{
  gid[i] <- gene.names[as.character(go.info$X.proteinId[i])]
  go[i] <- go.info$goAcc[i]
}

go.data.frame <- data.frame(GID=gid,GO=go,EVIDENCE=evidence)
head(go.data.frame)

go.data.frame <- data.frame(GID=go.info$X.proteinId,GO=go.info$goAcc,EVIDENCE=evidence)

length(go.data.frame$GID)
length(unique(go.data.frame$GID))

## For the rest of the data frames. According to the help info "the 1st column of every 
## data.frame must be labeled GID, and correspond to a gene ID that is universal for the 
## entire set of data.frames.  The GID is how the different tables will be joined internally"

## Generate enzyme data.frame

kegg.info <- read.table(file = "Ostta4221_3_GeneCatalog_proteins_20161028_KEGG.tab",header = T,sep = "\t",as.is = T,comment.char = "")
head(kegg.info)

gid <- vector(mode="character",length=nrow(kegg.info))
enzyme <- vector(mode="character",length=nrow(kegg.info))

for(i in 1:nrow(kegg.info))
{
  gid[i] <- gene.names[as.character(kegg.info$X.proteinId[i])]
  enzyme[i] <- kegg.info$ecNum[i]
}

enzyme.data.frame <- data.frame(GID=gid,ENZYME=enzyme)
head(enzyme.data.frame)

enzyme.data.frame <- data.frame(GID=kegg.info$X.proteinId,ENZYME=kegg.info$ecNum)
head(enzyme.data.frame)


## Generate KOG data.frame
kog.info <- read.table(file = "Ostta4221_3_GeneCatalog_proteins_20161028_KOG.tab",header = T,sep = "\t",as.is = T,comment.char = "")
head(kog.info)

gid <- vector(mode="character",length=nrow(kog.info))
kog <- vector(mode="character",length=nrow(kog.info))

for(i in 1:nrow(kog.info))
{
  gid[i] <- gene.names[as.character(kog.info$proteinId[i])]
  kog[i] <- kog.info$kogid[i]
}

kog.data.frame <- data.frame(GID=gid,KOG=kog)
head(kog.data.frame)

kog.data.frame <- data.frame(GID=kog.info$proteinId,KOG=kog.info$kogid)


## Ostreococcus tauri Taxonomy ID: 70448

## Load require package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationForge", version = "3.8")
library(AnnotationForge)

makeOrgPackage(go=go.data.frame,
#               ENZYME=enzyme.data.frame,
#               KOG=kog.data.frame,
#               SYMBOL=id.ostta.name,
               version = "0.1",
               maintainer = "Francisco J. Romero-Campero <fran@us.es>",
               author = "Francisco J. Romero-Campero",
               outputDir = ".", 
               tax_id = "70448",
               genus = "Ostreococcus",
               species = "tauri",
               goTable = "go",
               verbose = TRUE)

install.packages("./org.Otauri.eg.db/", repos=NULL)

library(org.Otauri.eg.db)
columns(org.Otauri.eg.db)
traceback()


data <- go.data.frame
unique(unlist(unname(lapply(data, "[", "GID"))))
unlist(unname(lapply(data, "[", "GID")))
structure(res, levels = lv, names = nm, class = "factor")

library(org.At.tair.db)
columns(org.At.tair.db)
head(select(org.At.tair.db,columns = c("GO"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("ENZYME"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("PATH"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("GENENAME"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("SYMBOL"),keys=keys(org.At.tair.db,keytype = "TAIR")))


head(select(org.At.tair.db,columns = c("ENTREZID"),keys=keys(org.At.tair.db,keytype = "TAIR")))

head(go.data.frame)



######################################
######################################

## Makes an organism package for Zebra Finch data.frames:
finchFile <- system.file("extdata","finch_info.txt",package="AnnotationForge")
finch <- read.table(finchFile,sep="\t")

## not that this is how it should always be, but that it *could* be this way.
fSym <- finch[,c(2,3,9)]
fSym <- fSym[fSym[,2]!="-",]
fSym <- fSym[fSym[,3]!="-",]
colnames(fSym) <- c("GID","SYMBOL","GENENAME")

fChr <- finch[,c(2,7)]
fChr <- fChr[fChr[,2]!="-",]
colnames(fChr) <- c("GID","CHROMOSOME")

finchGOFile <- system.file("extdata","GO_finch.txt",package="AnnotationForge")
fGO <- read.table(finchGOFile,sep="\t")
fGO <- fGO[fGO[,2]!="",]
fGO <- fGO[fGO[,3]!="",]
colnames(fGO) <- c("GID","GO","EVIDENCE")
head(fGO)
typeof(fGO$GID)
head(go.data.frame)
makeOrgPackage(gene_info=fSym, chromosome=fChr, go=fGO,
               version="0.1",
               maintainer="Some One <so@someplace.org>",
               author="Some One <so@someplace.org>",
               outputDir = ".",
               tax_id="59729",
               genus="Taeniopygia",
               species="guttata",
               goTable="go")
