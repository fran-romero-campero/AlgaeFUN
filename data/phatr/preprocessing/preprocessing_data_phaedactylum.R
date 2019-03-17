## Author: Francisco J. Romero-Campero, Ana Belén Romero-Losada
## Date: March 2019
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

## Read fasta file and split sequence names with token | for transcripts
transcripts.info <- read.fasta(file="Ostta4221_3_GeneCatalog_transcripts_20161028.nt.fasta",seqtype = "DNA")
trans.names <- getName(transcripts.info)
splitted.trans.names <- strsplit(trans.names,split="\\|")

## Loop to extract protein ids and ostta names. Initialize accumulators
transcripts.id <- vector(mode="character",length=length(trans.names))
ostta.trans.names <- vector(mode="character",length=length(trans.names))

for(i in 1:length(trans.names))
{
  transcripts.id[i] <- splitted.trans.names[[i]][3]
  ostta.trans.names[i] <- substr(x = splitted.trans.names[[i]][4],start = 4,stop = 18)
}

## Generate and write output data frame
id.ostta.name <- data.frame(GID=proteins.id,GENENAME=ostta.names,stringsAsFactors = FALSE)
head(id.ostta.name)

write.table(x = id.ostta.name,file = "correspondence_gene_id_ostta_name.tsv",sep="\t",quote=F,row.names = F)

id.trans.ostta.name <- data.frame(TID=transcripts.id,TRANSNAME=ostta.trans.names,stringsAsFactors = F)
write.table(x = id.trans.ostta.name, file="correspondence_trans_id_ostta_name.tsv",sep="\t",quote=F,row.names = F)

id.ostta.name <- read.table(file = "correspondence_gene_id_ostta_name.tsv",sep="\t",header=T,as.is=T)
head(id.ostta.name)

gene.ids <- id.ostta.name$GID
gene.names <- id.ostta.name$GENENAME
names(gene.names) <- gene.ids

## Generate SYMBOL data frame
symbol.data.frame <- data.frame(GID=gene.names,SYMBOL=gene.names,stringsAsFactors = FALSE)

## Remove duplicated rows
symbol.data.frame <- symbol.data.frame[!duplicated(symbol.data.frame),]

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


go.data.frame <- data.frame(GID=gid,GO=go,EVIDENCE=evidence,stringsAsFactors = FALSE)
head(go.data.frame)

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
 
enzyme.data.frame <- data.frame(GID=gid,ENZYME=enzyme,stringsAsFactors = FALSE)
head(enzyme.data.frame)

## Remove duplicated rows
enzyme.data.frame <- enzyme.data.frame[!duplicated(enzyme.data.frame),]
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

kog.data.frame <- data.frame(GID=gid,KOG=kog,stringsAsFactors = FALSE)
head(kog.data.frame)

## Remove duplicated rows
kog.data.frame <- kog.data.frame[!duplicated(kog.data.frame),]

## Ostreococcus tauri Taxonomy ID: 70448

## Load require package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationForge", version = "3.8")
library(AnnotationForge)

makeOrgPackage(go=go.data.frame,
               SYMBOL=symbol.data.frame,
               ENZYME=enzyme.data.frame,
               KOG=kog.data.frame,
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
head(select(org.Otauri.eg.db,columns = c("GO"),keys=keys(org.Otauri.eg.db,keytype = "GID")))
head(select(org.Otauri.eg.db,columns = c("ENZYME"),keys=keys(org.Otauri.eg.db,keytype = "GID")))
head(select(org.Otauri.eg.db,columns = c("KOG"),keys=keys(org.Otauri.eg.db,keytype = "GID")))

library(org.At.tair.db)
columns(org.At.tair.db)
head(select(org.At.tair.db,columns = c("GO"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("ENZYME"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("PATH"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("GENENAME"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("SYMBOL"),keys=keys(org.At.tair.db,keytype = "TAIR")))


ostta.example <- read.table(file = "../ota_trough_dark_light_peak_light_dark.txt",header = FALSE,as.is = TRUE)[[1]]
length(ostta.example)
BiocManager::install("clusterProfiler", version = "3.8")
library(clusterProfiler)

ostta.universe <- unique(select(org.Otauri.eg.db,columns = c("GO"),keys=keys(org.Otauri.eg.db,keytype = "GID"))[["GID"]])
length(ostta.universe)

keytypes(org.Otauri.eg.db)

ego <- enrichGO(gene          = ostta.example,
                universe      = ostta.universe,
                OrgDb         = org.Otauri.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE,
                keyType = "GID")

barplot(ego,drop=TRUE,showCategory = 10)
goplot(ego)
dotplot(ego)
emapplot(ego)
cnetplot(ego)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview", version = "3.8")
library("pathview")

help(enrichKEGG)

kk <- enrichKEGG(gene = paste0("OT_",ostta.example), organism = "ota",keyType = "kegg",
                 universe = paste0("OT_",ostta.universe),qvalueCutoff = 0.05)

head(kk)


mkk <- enrichMKEGG(gene = paste0("OT_",ostta.example), organism = "ota",keyType = "kegg")
head(mkk)

genes.pathway <- rep(0,length(ostta.universe))
names(genes.pathway) <- paste0("OT_",ostta.universe)

genes.pathway[paste0("OT_",ostta.example)] <- 1

pathview(gene.data = sort(genes.pathway,decreasing = TRUE), pathway.id = "ota03030", species = "ota",limit = list(gene=max(abs(genes.pathway)), cpd=1),gene.idtype ="kegg")

## Preprocess gff3
library(stringr)
ostta.gff3.0 <- read.table(file="ostreococcus_tauri.gff3",header=F,quote = "#",as.is=T)
ostta.gff3 <- ostta.gff3.0
head(ostta.gff3)


# Correspondence between gene/transcript names and ids
gene.names <- id.ostta.name$GENENAME
names(gene.names) <- id.ostta.name$GID

trans.names <- id.trans.ostta.name$TRANSNAME
names(trans.names) <- id.trans.ostta.name$TID

## Replace ids with names
for(i in 1:nrow(ostta.gff3))
{
  current.element <- ostta.gff3$V9[i]
  
  current.protein.id  <- strsplit(strsplit(current.element,split="proteinId=")[[1]][2],split=";")[[1]][1]
  
  # Some rows do not represent a transcript and do not need any replacement
  if(!is.na(current.protein.id))
  {
    current.protein.name <- gene.names[[current.protein.id]]
    
    current.trans.id <- strsplit(strsplit(current.element,split="transcriptId=")[[1]][2],split=";")[[1]][1]
    current.trans.name <- trans.names[[current.trans.id]]
    
    current.element <- str_replace_all(string = current.element,pattern = current.protein.id,replacement = current.protein.name)
    current.element <- str_replace_all(string = current.element,pattern = current.trans.id,replacement = current.trans.name)
    
    ostta.gff3$V9[i] <- current.element
  }
}

write.table(x = ostta.gff3,file = "ostreococcus_tauri_annotation.gff3",quote = F,sep = "\t",row.names = F,col.names = F)

## Manually edit the previous file to add ##gff-version 3
## Manually remove jgi.p|Ostta4221_3|

## Generate TxDb package from gff3 file
library("GenomicFeatures")

## Generate chromosome.info data frame
library(seqinr)

ostta.genome.data <- read.fasta(file = "ostreococcus_tauri_genome.fasta",seqtype = "DNA")
chromosome.names <- getName(ostta.genome.data)
chromosome.lengths <- sapply(X=getSequence(ostta.genome.data),FUN = length)

chromosomes.info <- data.frame(chrom=chromosome.names,length=chromosome.lengths,is_circular=FALSE)

## Meta data info
meta.data.info <- data.frame(name=c("Resource URL","Genome"),value=c("https://genome.jgi.doe.gov/Ostta4221_3/Ostta4221_3.home.html","v3.0"))

ostta.txdb <- makeTxDbFromGFF(file = "ostreococcus_tauri_annotation.gff3",format = "gff3",dataSource = "JGI",organism = "Ostreococcus tauri",taxonomyId = 70448,chrominfo = chromosomes.info,metadata = meta.data.info)

ostta.txdb
genes(ostta.txdb)

makeTxDbPackage(txdb = ostta.txdb, version = "0.1", maintainer = "Francisco J. Romero-Campero <fran@us.es>", author = "Francisco J. Romero-Campero")

install.packages("./TxDb.Otauri.JGI/", repos=NULL)
## loading packages
#library(ChIPseeker)
library(TxDb.Otauri.JGI)
txdb <- TxDb.Otauri.JGI
genes(txdb)
