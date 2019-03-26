## Author: Francisco J. Romero-Campero, Ana Belén Romero-Losada
## Date: February 2019
## Contact: fran@us.es

## Documentation on how to generate org.db packages can be found on these links:
## http://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html
## 
vocar.info <- read.table(file="Vcarteri_317_v2.1.annotation_info.txt",header=T,comment.char = "",sep = "\t",as.is=T, fill=T)
head(vocar.info)
nrow(vocar.info)

write.table(x = unique(vocar.info$locusName),file = "vocar_universe.txt",quote = F,row.names = F,col.names = F)



## Generate and write output data frame
id.vocar.name <- data.frame(GID=vocar.info$locusName,GENENAME=vocar.info$locusName,stringsAsFactors = FALSE)
head(id.vocar.name)

"Vocar.0012s0145" %in% id.vocar.name$GID

## Generate SYMBOL data frame
symbol.data.frame <- data.frame(GID=vocar.info$locusName,
                                SYMBOL=vocar.info$locusName,stringsAsFactors = FALSE)
head(symbol.data.frame)
nrow(symbol.data.frame)
symbol.data.frame[1:20,]

## Remove duplicated rows
symbol.data.frame <- symbol.data.frame[!duplicated(symbol.data.frame),]
nrow(symbol.data.frame)
sum(duplicated(symbol.data.frame))
head(symbol.data.frame)

"Vocar.0012s0145" %in% symbol.data.frame$GID


## Generate GO data.frame to create org.Db package
## Before reading the GO annotation file I have to replace the symbol ' by the word prime
## to avoid the error "EOF within quoted string"

## According to the vignettes of AnnotationForge bioconductor package for "Making Organisms
## Packages: "However to use the goTable argument, you have to follow a strict convention with 
## the data. Such a data.frame must have three columns only and these must correspond to the 
## gene id, GO id and evidence codes. These columns also have to be named as “GID”, “GO” and “EVIDENCE”

gid <- vocar.info$locusName
go <- vocar.info$GO
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
go.data.frame[1:10,]
nrow(go.data.frame)
length(unique(go.data.frame$GID))
sum(duplicated(go.data.frame))

"Vocar.0012s0145" %in% go.data.frame$GID

## For the rest of the data frames. According to the help info "the 1st column of every 
## data.frame must be labeled GID, and correspond to a gene ID that is universal for the 
## entire set of data.frames.  The GID is how the different tables will be joined internally"

## Generate enzyme data.frame

enzyme <- vocar.info$KEGG.ec
 
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

"Vocar.0012s0145" %in% enzyme.data.frame$GID
subset(enzyme.data.frame, GID == "Vocar.0012s0145")

## Generate KOG data.frame
kog <- vocar.info$KOG

# for(i in 1:nrow(kog.info))
# {
#   gid[i] <- gene.names[as.character(kog.info$proteinId[i])]
#   kog[i] <- kog.info$kogid[i]
# }

kog.data.frame <- data.frame(GID=gid,KOG=kog,stringsAsFactors = FALSE)
head(kog.data.frame)

## Remove duplicated rows
kog.data.frame <- kog.data.frame[!duplicated(kog.data.frame),]
sum(duplicated(kog.data.frame))
nrow(kog.data.frame)
head(kog.data.frame)
length(unique(kog.data.frame$GID))


"Vocar.0012s0145" %in% kog.data.frame$GID

##Generate KO data frame
ko <- vocar.info$KO

ko.data.frame <- data.frame(GID=gid, KO=ko, stringsAsFactors = F)
head(ko.data.frame)
ko.data.frame[1:20,]

ko.data.frame <- ko.data.frame[!duplicated(ko.data.frame),]
nrow(ko.data.frame)
sum(duplicated(ko.data.frame))

"Vocar.0012s0145" %in% ko.data.frame$GID


##Generate panther data frame

panther <- vocar.info$Panther
panther.data.frame <- data.frame(GID=gid, PANTHER=panther, stringsAsFactors = F)
nrow(panther.data.frame)
head(panther.data.frame)

panther.data.frame <- panther.data.frame[!duplicated(panther.data.frame),]
nrow(panther.data.frame)
sum(duplicated(panther.data.frame))
#panther.expanded.data.frame <- c()
# for(i in 1:nrow(panther.data.frame))
# {
#   panther.splitted <- strsplit(panther.data.frame$PANTHER[i],split=",")[[1]]
#   
#   if(length(panther.splitted) > 0)
#   {
#     for(j in 1:length(panther.splitted))
#     {
#       panther.expanded.data.frame <- rbind(panther.expanded.data.frame,c(panther.data.frame[i,1],panther.splitted[j]))
#     }
#   } else
#   {
#     panther.expanded.data.frame <- rbind(panther.expanded.data.frame,c(panther.data.frame[i,1],""))
#   }
# }
# head (panther.expanded.data.frame)
panther.data.frame <- data.frame(GID=panther.data.frame[,1],
                           PANTHER=panther.data.frame[,2],
                          stringsAsFactors = FALSE)
head(panther.data.frame)

"Vocar.0012s0145" %in% panther.data.frame$GID

## Dunaliella salina Taxonomy ID: 3046

## Load require package
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("AnnotationForge", version = "3.8")
# BiocManager::install("GO.db", version = "3.8")
library(AnnotationForge)
library(GO.db)

"Vocar.0012s0145" %in% go.data.frame$GID
"Vocar.0012s0145" %in% symbol.data.frame$GID
"Vocar.0012s0145" %in% enzyme.data.frame$GID
"Vocar.0012s0145" %in% kog.data.frame$GID
"Vocar.0012s0145" %in% ko.data.frame$GID
"Vocar.0012s0145" %in% panther.data.frame$GID

makeOrgPackage(go=go.data.frame,
               SYMBOL=symbol.data.frame,
               ENZYME=enzyme.data.frame,
               KOG=kog.data.frame,
               KO=ko.data.frame,
               PANTHER=panther.data.frame,
               version = "0.1",
               maintainer = "Francisco J. Romero-Campero <fran@us.es>",
               author = "Ana B. Romero-Losada",
               outputDir = ".", 
               tax_id = "3067",
               genus = "Volvox",
               species = "carteri",
               goTable = "go",
               verbose = TRUE)

install.packages("./org.Vcarteri.eg.db/", repos=NULL)
remove.packages("org.Vcarteri.eg.db")

library(org.Vcarteri.eg.db)
columns(org.Vcarteri.eg.db)
head(select(org.Vcarteri.eg.db,columns = c("ENZYME"),keys=keys(org.Vcarteri.eg.db,keytype = "GID")))
head(select(org.Vcarteri.eg.db,columns = c("EVIDENCE"),keys=keys(org.Vcarteri.eg.db,keytype = "GID")))
head(select(org.Vcarteri.eg.db,columns = c("EVIDENCEALL"),keys=keys(org.Vcarteri.eg.db,keytype = "GID")))
head(select(org.Vcarteri.eg.db,columns = c("GO"),keys=keys(org.Vcarteri.eg.db,keytype = "GID")))
head(select(org.Vcarteri.eg.db,columns = c("GOALL"),keys=keys(org.Vcarteri.eg.db,keytype = "GID")))
head(select(org.Vcarteri.eg.db,columns = c("KOG"),keys=keys(org.Vcarteri.eg.db,keytype = "GID")))
head(select(org.Vcarteri.eg.db,columns = c("ONTOLOGY"),keys=keys(org.Vcarteri.eg.db,keytype = "GID")))
head(select(org.Vcarteri.eg.db,columns = c("ONTOLOGYALL"),keys=keys(org.Vcarteri.eg.db,keytype = "GID")))
head(select(org.Vcarteri.eg.db,columns = c("SYMBOL"),keys=keys(org.Vcarteri.eg.db,keytype = "GID")))
head(select(org.Vcarteri.eg.db,columns = c("KO"),keys=keys(org.Vcarteri.eg.db,keytype = "GID")))
head(select(org.Vcarteri.eg.db,columns = c("PANTHER"),keys=keys(org.Vcarteri.eg.db,keytype = "GID")))


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
vocar.gff3 <- read.table(file="../annotation/volvox_carteri.gff3",header=F,quote = "#",as.is=T)
vocar.gtf <- vocar.gff3
head(vocar.gff3)

unique(vocar.gff3$V3)

extract.attributes <- function(str.attributes)
{
  split.attribs <- strsplit(str.attributes,split=";")[[1]]
  
  attrib.names <- vector(mode="character", length=length(split.attribs))
  attrib.values <- vector(mode="character", length=length(split.attribs))

  for(i in 1:length(split.attribs))
  {
    attrib.names[i] <-strsplit(split.attribs[i],split="=")[[1]][1]
    attrib.values[i] <- strsplit(split.attribs[i],split="=")[[1]][2]
  }
  
  names(attrib.values) <- attrib.names
  
  return(attrib.values)
}

for(i in 1:nrow(vocar.gff3))
{
  current.attributes <- extract.attributes(vocar.gff3$V9[i])
  if(vocar.gff3$V3[i] == "gene")
  {
    gene.id <- current.attributes["Name"]
    vocar.gtf$V9[i] <- paste("gene_id", paste("\"",gene.id,"\";",sep=""))
  } else if(vocar.gff3$V3[i] == "mRNA")
  {
    gene.id <- substr(x = current.attributes["Name"], start = 1, stop = 15)
    transcript.id <- current.attributes["Name"]
    vocar.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep="")), collapse = " ")
  } else if(vocar.gff3$V3[i] == "exon")
  {
    gene.id <- substr(x = current.attributes["ID"], start = 1, stop = 15) #substr(strsplit(current.attributes[1],split="=")[[1]][2],start = 1,stop = 13)
    trancript.id <- substr(x = current.attributes["ID"], start = 1, stop = 17) #substr(strsplit(current.attributes[1],split="=")[[1]][2],start = 1,stop = 18)
    exon.number <- strsplit(current.attributes["ID"],split="exon.")[[1]][2]
    vocar.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep=""),"exon_number",paste("\"",exon.number,"\";",sep="")), collapse = " ")
  } else if(vocar.gff3$V3[i] == "five_prime_UTR" || vocar.gff3$V3[i] == "three_prime_UTR" || vocar.gff3$V3[i] == "CDS")
  {
    gene.id <- substr(x = current.attributes["ID"], start = 1, stop = 15) # substr(strsplit(current.attributes[2],split="=")[[1]][2],start = 1,stop = 13)
    transcript.id <- substr(x = current.attributes["ID"], start = 1, stop = 17) #substr(strsplit(current.attributes[2],split="=")[[1]][2],start = 1,stop = 18)
    vocar.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep="")), collapse = " ")
  }
}

head(vocar.gtf)

write.table(x = vocar.gtf,file = "volvox_carteri.gtf",sep = "\t",row.names = F,col.names = F,quote = F)

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


## Processing of fasta files with proteins
library(seqinr)

proteins.data <- read.fasta(file="GCF_000143455_1_v1_0_translated_cds.fa",seqtype = "AA")

proteins.annot <- getAnnot(proteins.data)
proteins.seqs <- getSequence(proteins.data)

extract.volcadraft <- function(gene.annotation)
{
  return(strsplit(strsplit(gene.annotation,split="locus_tag=")[[1]][2],split="]")[[1]][1])  
}

new.names <- sapply(proteins.annot,extract.volcadraft)
is.vector(new.names)

write.fasta(sequences = proteins.seqs, names = new.names, file.out = "volcadraft_seqs.fasta")

volca.data <- read.fasta(file = "Vcarteri_317_v2_1_protein.fa",seqtype = "AA")

volca.names <- getName(volca.data)

extract.volca <- function(gene.name)
{
  return(substr(x = gene.name,start = 1,stop = 15))  
}

volca.new.names <- sapply(X = volca.names, FUN = extract.volca)
names(volca.new.names) <- NULL

write.fasta(sequences = getSequence(volca.data), names = volca.new.names,file.out = "vcarteri_protein.fa")
