## Author: Francisco J. Romero-Campero, Ana Belén Romero-Losada
## Date: March 2019
## Contact: fran@us.es

## Documentation on how to generate org.db packages can be found on these links:
## http://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html
## 

## Load correspondence between naga ids from ensembl and jgi
naga.ensembl.jgi <- read.table(file="jgi_vs_ensembl/result_identity.txt",as.is=T)
head(naga.ensembl.jgi)

naga.ensembl <- naga.ensembl.jgi$V1
naga.jgi <- naga.ensembl.jgi$V2

names(naga.jgi) <- naga.ensembl
names(naga.ensembl) <- naga.jgi

gene.names <- unique(naga.ensembl)
sum(is.na(gene.names))
write.table(x = gene.names, file = "naga_universe.txt",quote = F, row.names = F, col.names = F)

naga.go.info <- read.table(file="Nangad1_GeneCatalog_proteins_20180307_GO.tab",header=T,sep="\t",as.is=T)
nrow(naga.go.info)
head(naga.go.info)

## Generate SYMBOL data frame
symbol.data.frame <- data.frame(GID=gene.names,SYMBOL=gene.names,stringsAsFactors = FALSE)
nrow(symbol.data.frame)
## Remove duplicated rows
symbol.data.frame <- symbol.data.frame[!duplicated(symbol.data.frame),]
nrow(symbol.data.frame)
head(symbol.data.frame)

## Generate GO data.frame to create org.Db package
## Before reading the GO annotation file I have to replace the symbol ' by the word prime
## to avoid the error "EOF within quoted string"

## According to the vignettes of AnnotationForge bioconductor package for "Making Organisms
## Packages: "However to use the goTable argument, you have to follow a strict convention with 
## the data. Such a data.frame must have three columns only and these must correspond to the 
## gene id, GO id and evidence codes. These columns also have to be named as “GID”, “GO” and “EVIDENCE”

gid <- vector(mode="character",length=nrow(naga.go.info))
go <- vector(mode="character",length=nrow(naga.go.info))
evidence <- rep("ISS",nrow(naga.go.info))

for(i in 1:nrow(naga.go.info))
{
  gid[i] <- naga.ensembl[as.character(naga.go.info$proteinId[i])]
  go[i] <- naga.go.info$goAcc[i]
}

go.data.frame <- data.frame(GID=gid,GO=go,EVIDENCE=evidence,stringsAsFactors = FALSE)
go.data.frame <- go.data.frame[!duplicated(go.data.frame),]
go.data.frame <- go.data.frame[go.data.frame$GO != "",]
go.data.frame <- go.data.frame[!is.na(go.data.frame$GID),]
head(go.data.frame)

## For the rest of the data frames. According to the help info "the 1st column of every 
## data.frame must be labeled GID, and correspond to a gene ID that is universal for the 
## entire set of data.frames.  The GID is how the different tables will be joined internally"

## Generate enzyme data.frame
naga.ec.info <- read.table(file="Nangad1_GeneCatalog_proteins_20180307_KEGG.tab",header=T,sep="\t",as.is=T)
nrow(naga.ec.info)
head(naga.ec.info)

gid <- vector(mode="character",length=nrow(naga.ec.info))
enzyme <- vector(mode="character",length=nrow(naga.ec.info))

for(i in 1:nrow(naga.ec.info))
{
  gid[i] <- naga.ensembl[as.character(naga.ec.info$proteinId[i])] 
  enzyme[i] <- naga.ec.info$ecNum[i]
}

enzyme.data.frame <- data.frame(GID=gid,ENZYME=enzyme,stringsAsFactors = FALSE)
head(enzyme.data.frame)

## Remove duplicated rows
enzyme.data.frame <- enzyme.data.frame[!duplicated(enzyme.data.frame),]
enzyme.data.frame <- enzyme.data.frame[!is.na(enzyme.data.frame$ENZYME),]
enzyme.data.frame <- enzyme.data.frame[!is.na(enzyme.data.frame$GID),]
head(enzyme.data.frame)
nrow(enzyme.data.frame)
enzyme.data.frame$ENZYME

## Generate KOG data.frame
naga.kog.info <- read.table(file="Nangad1_GeneCatalog_proteins_20180307_KOG.tab",header=T,sep="\t",as.is=T)
nrow(naga.kog.info)
head(naga.kog.info)

gid <- vector(mode="character",length=nrow(naga.kog.info))
kog <- vector(mode="character",length=nrow(naga.kog.info))

for(i in 1:nrow(naga.kog.info))
{
  gid[i] <- naga.ensembl[as.character(naga.kog.info$proteinId[i])] 
  kog[i] <- naga.kog.info$kogid[i]
}

kog.data.frame <- data.frame(GID=gid,KOG=kog,stringsAsFactors = FALSE)
head(kog.data.frame)

## Remove duplicated rows
kog.data.frame <- kog.data.frame[!duplicated(kog.data.frame),]
kog.data.frame <- kog.data.frame[!is.na(kog.data.frame$KOG),]
kog.data.frame <- kog.data.frame[!is.na(kog.data.frame$GID),]
head(kog.data.frame)
nrow(kog.data.frame)


## Naga draft for kegg analysis
naga.draft <- read.table(file="ensembl_vs_kegg/naga_nga_correspondence.txt",header=F,as.is=T)
head(naga.draft)
naga.draft$V1
naga.draft$V2

nagadraft.data.frame <- data.frame(GID=naga.draft$V1,NAGADRAFT=naga.draft$V2,stringsAsFactors = FALSE)
head(nagadraft.data.frame)
nagadraft.data.frame <- nagadraft.data.frame[!duplicated(nagadraft.data.frame),]
nrow(nagadraft.data.frame)

## Nannochloropsis gaditana Taxonomy ID: 1093141

## Load require package
library(AnnotationForge)
naga.go.data.frame <- go.data.frame
naga.go.data.frame$GID

head(naga.go.data.frame)
head(symbol.data.frame)
head(enzyme.data.frame)
head(kog.data.frame)
head(nagadraft.data.frame)

makeOrgPackage(go=go.data.frame,
               SYMBOL=symbol.data.frame,
               ENZYME=enzyme.data.frame,
               KOG=kog.data.frame,
               NAGADRAFT=nagadraft.data.frame,
               version = "0.1",
               maintainer = "Francisco J. Romero-Campero <fran@us.es>",
               author = "Francisco J. Romero-Campero",
               outputDir = ".", 
               tax_id = "1093141",
               genus = "Nannochloropsis",
               species = "gaditana",
               goTable = "go",
               verbose = TRUE)

install.packages("./org.Ngaditana.eg.db/", repos=NULL)
#remove.packages("org.Ptricornutum.eg.db")
library(org.Ngaditana.eg.db)
columns(org.Ngaditana.eg.db)
head(select(org.Ngaditana.eg.db,columns = c("GO"),keys=keys(org.Ngaditana.eg.db,keytype = "GID")))
head(select(org.Ngaditana.eg.db,columns = c("SYMBOL"),keys=keys(org.Ngaditana.eg.db,keytype = "GID")))
head(select(org.Ngaditana.eg.db,columns = c("ENZYME"),keys=keys(org.Ngaditana.eg.db,keytype = "GID")))
head(select(org.Ngaditana.eg.db,columns = c("NAGADRAFT"),keys=keys(org.Ngaditana.eg.db,keytype = "GID")))

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


ptricornutum.gtf <- read.table(file = "../annotation/phaeodactylum_tricornutum.gtf",header = F,sep = "\t",as.is=T)
head(ptricornutum.gtf)

i <- 3

genes.in.gtf <- vector(mode="character", length=nrow(ptricornutum.gtf))
for(i in 1:nrow(ptricornutum.gtf))
{
  attrib.gtf <- ptricornutum.gtf$V9[i]
  genes.in.gtf[i] <- strsplit(strsplit(x = attrib.gtf, split="gene_id")[[1]][2], split=";")[[1]][1]
}

genes.in.gtf <- unique(genes.in.gtf)
length(genes.in.gtf)
genes.in.universe <- read.table(file = "../../../web_app/universe/phatri_universe.txt",as.is=T,header = F)[[1]]
length(genes.in.universe)

gene.universe.final <- unique(c(genes.in.gtf, genes.in.universe))
length(c(genes.in.gtf, genes.in.universe))
length(gene.universe.final)

intersect(genes.in.universe, genes.in.gtf)
genes.in.universe[2]
genes.in.gtf[2]

"Phatr3_J43104" %in% genes.in.universe
"Phatr3_J43104" %in% genes.in.gtf
genes.in.gtf[7397]


remove.white.spaces <- function(x)
{
  return(gsub(" ", "", x, fixed = TRUE))
}

genes.in.gtf <- sapply(genes.in.gtf,FUN = remove.white.spaces)
names(genes.in.gtf) <- NULL
genes.in.gtf

gene.universe.final <- unique(c(genes.in.gtf, genes.in.universe))
length(gene.universe.final)

write.table(x = gene.universe.final,file = "../../../web_app/universe/phatri_universe.txt",row.names = F, col.names = F)
