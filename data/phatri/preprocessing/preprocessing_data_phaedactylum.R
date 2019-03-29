## Author: Francisco J. Romero-Campero, Ana Belén Romero-Losada
## Date: March 2019
## Contact: fran@us.es

## Documentation on how to generate org.db packages can be found on these links:
## http://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html
## 

phatri.info <- read.table(file="phatri_info.txt",header=T,sep="\t",as.is=T)
nrow(phatri.info)
gene.names <- unique(phatri.info$Gene.stable.ID)
write.table(x = gene.names, file = "phatri_universe.txt",quote = F, row.names = F, col.names = F)

## Generate SYMBOL data frame
symbol.data.frame <- data.frame(GID=gene.names,SYMBOL=gene.names,stringsAsFactors = FALSE)
nrow(symbol.data.frame)
## Remove duplicated rows
symbol.data.frame <- symbol.data.frame[!duplicated(symbol.data.frame),]
nrow(symbol.data.frame)

## Generate GO data.frame to create org.Db package
## Before reading the GO annotation file I have to replace the symbol ' by the word prime
## to avoid the error "EOF within quoted string"

## According to the vignettes of AnnotationForge bioconductor package for "Making Organisms
## Packages: "However to use the goTable argument, you have to follow a strict convention with 
## the data. Such a data.frame must have three columns only and these must correspond to the 
## gene id, GO id and evidence codes. These columns also have to be named as “GID”, “GO” and “EVIDENCE”

gid <- vector(mode="character",length=nrow(phatri.info))
go <- vector(mode="character",length=nrow(phatri.info))
evidence <- rep("ISS",nrow(phatri.info))

for(i in 1:nrow(phatri.info))
{
  gid[i] <- phatri.info$Gene.stable.ID[i]   #gene.names[as.character(go.info$X.proteinId[i])]
  go[i] <- phatri.info$GO.term.accession[i]
}


go.data.frame <- data.frame(GID=gid,GO=go,EVIDENCE=evidence,stringsAsFactors = FALSE)
go.data.frame <- go.data.frame[!duplicated(go.data.frame),]
go.data.frame <- go.data.frame[go.data.frame$GO != "",]
head(go.data.frame)


## For the rest of the data frames. According to the help info "the 1st column of every 
## data.frame must be labeled GID, and correspond to a gene ID that is universal for the 
## entire set of data.frames.  The GID is how the different tables will be joined internally"

## Generate enzyme data.frame

gid <- vector(mode="character",length=nrow(phatri.info))
enzyme <- vector(mode="character",length=nrow(phatri.info))

for(i in 1:nrow(phatri.info))
{
  gid[i] <- phatri.info$Gene.stable.ID[i] 
  enzyme[i] <- strsplit(phatri.info$KEGG.Pathway.and.Enzyme.ID[i], split="\\+")[[1]][2]
}

enzyme.data.frame <- data.frame(GID=gid,ENZYME=enzyme,stringsAsFactors = FALSE)
head(enzyme.data.frame)

## Remove duplicated rows
enzyme.data.frame <- enzyme.data.frame[!duplicated(enzyme.data.frame),]
enzyme.data.frame <- enzyme.data.frame[!is.na(enzyme.data.frame$ENZYME),]
head(enzyme.data.frame)
enzyme.data.frame$ENZYME

##Generate panther data frame
panther <- phatri.info$Gene.stable.ID
panther.data.frame <- data.frame(GID=phatri.info$Gene.stable.ID, PANTHER=phatri.info$PANTHER.ID, stringsAsFactors = F)
nrow(panther.data.frame)
head(panther.data.frame)

panther.data.frame <- panther.data.frame[!duplicated(panther.data.frame),]
nrow(panther.data.frame)
panther.data.frame <- panther.data.frame[panther.data.frame$PANTHER != "",]
head(panther.data.frame)

##Generate pfam data frame
pfam.data.frame <- data.frame(GID=phatri.info$Gene.stable.ID, PFAM=phatri.info$Pfam.ID, stringsAsFactors = F)
nrow(pfam.data.frame)
head(pfam.data.frame)

pfam.data.frame <- pfam.data.frame[!duplicated(pfam.data.frame),]
nrow(pfam.data.frame)
pfam.data.frame <- pfam.data.frame[pfam.data.frame$PFAM != "",]
head(pfam.data.frame)

## Phaeodactylum tricornutum Taxonomy ID: 556484

## Load require package
library(AnnotationForge)

makeOrgPackage(go=go.data.frame,
               SYMBOL=symbol.data.frame,
               ENZYME=enzyme.data.frame,
               PANTHER=panther.data.frame,
               PFAM=pfam.data.frame,
               version = "0.1",
               maintainer = "Francisco J. Romero-Campero <fran@us.es>",
               author = "Francisco J. Romero-Campero",
               outputDir = ".", 
               tax_id = "556484",
               genus = "Phaeodactylum",
               species = "tricornutum",
               goTable = "go",
               verbose = TRUE)

install.packages("./org.Ptricornutum.eg.db", repos=NULL)
remove.packages("org.Ptricornutum.eg.db")
library(org.Ptricornutum.eg.db)
columns(org.Ptricornutum.eg.db)
head(select(org.Ptricornutum.eg.db,columns = c("GO"),keys=keys(org.Ptricornutum.eg.db,keytype = "GID")))
head(select(org.Ptricornutum.eg.db,columns = c("SYMBOL"),keys=keys(org.Ptricornutum.eg.db,keytype = "GID")))
head(select(org.Ptricornutum.eg.db,columns = c("ENZYME"),keys=keys(org.Ptricornutum.eg.db,keytype = "GID")))
head(select(org.Ptricornutum.eg.db,columns = c("PFAM"),keys=keys(org.Ptricornutum.eg.db,keytype = "GID")))

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
