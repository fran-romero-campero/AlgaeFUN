## Author: Francisco J. Romero-Campero, Ana Bel??n Romero-Losada
## Date: November 2018
## Contact: fran@us.es; ABRomeroLosada@gmail.com

## Documentation on how to generate org.db packages can be found on these links:
## http://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html
## 

salina.info <- read.table(file="Dsalina_325_v1.0.annotation_info.txt",header=T,comment.char = "",sep = "\t",as.is=T, quote = "")
head(salina.info)

## Generate and write output data frame
id.salina.name <- data.frame(GID=salina.info$locusName,GENENAME=salina.info$locusName,stringsAsFactors = FALSE)
head(id.salina.name)

## Generate SYMBOL data frame
symbol.data.frame <- data.frame(GID=salina.info$locusName,SYMBOL=salina.info$locusName,stringsAsFactors = FALSE)
head(symbol.data.frame)

## Remove duplicated rows
symbol.data.frame <- symbol.data.frame[!duplicated(symbol.data.frame),]

## Generate GO data.frame to create org.Db package
## Before reading the GO annotation file I have to replace the symbol ' by the word prime
## to avoid the error "EOF within quoted string"

## According to the vignettes of AnnotationForge bioconductor package for "Making Organisms
## Packages: "However to use the goTable argument, you have to follow a strict convention with 
## the data. Such a data.frame must have three columns only and these must correspond to the 
## gene id, GO id and evidence codes. These columns also have to be named as ???GID???, ???GO??? and ???EVIDENCE???

gid <- salina.info$locusName
go <- salina.info$GO
evidence <- rep("ISS",length(gid))

#Create the GO data frame
go.data.frame <- data.frame(GID=gid,GO=go,EVIDENCE=evidence,stringsAsFactors = FALSE)
head(go.data.frame)

#Locate each GO term in different rows with its respective gene id.

go.expanded.data.frame <- c()

for(i in 1:nrow(go.data.frame))
{
  go.splitted <- strsplit(go[i],split=",")[[1]]
  
  if(length(go.splitted) > 0)
  {
    for(j in 1:length(go.splitted))
    {
      go.expanded.data.frame <- rbind(go.expanded.data.frame,c(go.data.frame[i,1],go.splitted[j]))
    }
  }
}

head(go.expanded.data.frame)

go.data.frame <- data.frame(GID=go.expanded.data.frame[,1],
                            GO=go.expanded.data.frame[,2],
                            EVIDENCE=rep("ISS",nrow(go.expanded.data.frame)),stringsAsFactors = FALSE)
head(go.data.frame)


## Remove duplicated rows
go.data.frame <- go.data.frame[!duplicated(go.data.frame),]

## For the rest of the data frames. According to the help info "the 1st column of every 
## data.frame must be labeled GID, and correspond to a gene ID that is universal for the 
## entire set of data.frames.  The GID is how the different tables will be joined internally"

## Generate enzyme data.frame

enzyme <- salina.info$KEGG.ec

enzyme.data.frame <- data.frame(GID=gid,ENZYME=enzyme,stringsAsFactors = FALSE)
head(enzyme.data.frame)

## Remove duplicated rows
enzyme.data.frame <- enzyme.data.frame[!duplicated(enzyme.data.frame),]
head(enzyme.data.frame)

## Generate KOG data.frame
kog <- salina.info$KOG

kog.data.frame <- data.frame(GID=gid,KOG=kog,stringsAsFactors = FALSE)
head(kog.data.frame)

## Remove duplicated rows
kog.data.frame <- kog.data.frame[!duplicated(kog.data.frame),]

## Load require package
library(AnnotationForge)

makeOrgPackage(go=go.data.frame,
               SYMBOL=symbol.data.frame,
               ENZYME=enzyme.data.frame,
               KOG=kog.data.frame,
               version = "0.1",
               maintainer = "Francisco J. Romero-Campero, Ana Belen Romero-Losada",
               author = "Francisco J. Romero-Campero y Ana Belen Romero-Losada",
               outputDir = ".", 
               tax_id = "3046",
               genus = "Dunaliella",
               species = "salina",
               goTable = "go",
               verbose = TRUE)

install.packages("./org.Dsalina.eg.db/", repos=NULL)

library(org.Osalina.eg.db)
columns(org.Creinhardtii.eg.db)
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
