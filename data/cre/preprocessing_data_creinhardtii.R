## Author: Francisco J. Romero-Campero, Ana Belén Romero-Losada
## Date: November 2018
## Contact: fran@us.es

## Documentation on how to generate org.db packages can be found on these links:
## http://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html
## 

cre.info <- read.table(file="Creinhardtii_281_v5.6.annotation_info.txt",header=T,comment.char = "",sep = "\t",as.is=T)
head(cre.info)

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
i <- 1
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

library(org.Creinhardtii.eg.db)
columns(org.Creinhardtii.eg.db)
head(select(org.Creinhardtii.eg.db,columns = c("GO"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID")))
head(select(org.Creinhardtii.eg.db,columns = c("SYMBOL"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID")))

head(select(org.Otauri.eg.db,columns = c("ENZYME"),keys=keys(org.Otauri.eg.db,keytype = "GID")))
head(select(org.Otauri.eg.db,columns = c("KOG"),keys=keys(org.Otauri.eg.db,keytype = "GID")))

library(org.At.tair.db)
columns(org.At.tair.db)
head(select(org.At.tair.db,columns = c("GO"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("ENZYME"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("PATH"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("GENENAME"),keys=keys(org.At.tair.db,keytype = "TAIR")))
head(select(org.At.tair.db,columns = c("SYMBOL"),keys=keys(org.At.tair.db,keytype = "TAIR")))

all.go <- select(org.Creinhardtii.eg.db,columns = c("GO"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID"))
head(all.go)

starch.go <- "GO:0005975"

cre.example <- subset(all.go,GO == starch.go)$GID

cre.example <- read.table(file = "clusters/cre_trough_dark_light_peak_light_dark.txt",header = FALSE,as.is = TRUE)[[1]]

length(cre.example)

BiocManager::install("clusterProfiler", version = "3.8")
library(clusterProfiler)

cre.universe <- unique(select(org.Creinhardtii.eg.db,columns = c("GO"),keys=keys(org.Creinhardtii.eg.db,keytype = "GID"))[["GID"]])
length(cre.universe)

length(intersect(dna.example,cre.example))

keytypes(org.Creinhardtii.eg.db)

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
goplot(ego)
dotplot(ego)
emapplot(ego)
cnetplot(ego)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview", version = "3.8")
library("pathview")

help(enrichKEGG)

kk <- enrichKEGG(gene = cre.example, organism = "cre",keyType = "kegg",
                 universe = cre.universe,qvalueCutoff = 0.05)

draft.ids <- chlredraft.cre[cre.example]
names(draft.ids) <- NULL
draft.ids <- draft.ids[!is.na(draft.ids)]

draft.universe <- chlredraft.cre[cre.universe]
names(draft.universe) <- NULL
draft.universe <- draft.universe[!is.na(draft.universe)]

kk <- enrichKEGG(gene = draft.ids, organism = "cre",keyType = "kegg",
                 universe = draft.universe,qvalueCutoff = 0.05)




head(kk)


mkk <- enrichMKEGG(gene = draft.ids, organism = "cre",keyType = "kegg")
head(mkk)

genes.pathway <- rep(0,length(draft.universe))
names(genes.pathway) <- draft.universe

genes.pathway[draft.ids] <- 1

res.kk <- head(as.data.frame(kk))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE), pathway.id =res.kk$ID[4], species = "cre",limit = list(gene=max(abs(genes.pathway)), cpd=1),gene.idtype ="kegg")

dna.example <- c("Cre01.g003463","Cre01.g019750","Cre01.g029100","Cre01.g045850","Cre01.g017450","Cre01.g023150","Cre01.g025300","Cre01.g036050","Cre01.g053200","Cre01.g063632","Cre10.g424200","Cre10.g446400","Cre10.g453650","Cre10.g466550","Cre10.g423800","Cre10.g428433","Cre10.g446600","Cre10.g455850","Cre12.g556911","Cre12.g540927","Cre12.g553253","Cre12.g490150","Cre12.g512500","Cre12.g521200","Cre12.g524150","Cre12.g524350","Cre12.g534151","Cre13.g566900","Cre13.g604850","Cre13.g607500","Cre14.g619825","Cre16.g664301","Cre16.g670550","Cre16.g682950","Cre16.g685613","Cre17.g710150","Cre17.g718850","Cre17.g746347","Cre02.g079850","Cre02.g084800","Cre02.g082000","Cre03.g145687","Cre03.g158550","Cre03.g162250","Cre03.g175850","Cre03.g179550","Cre03.g181650","Cre03.g196600","Cre03.g202250","Cre03.g199400","Cre03.g172050","Cre03.g178650","Cre03.g179961","Cre03.g192550","Cre03.g204900","Cre04.g227000","Cre04.g214350","Cre04.g227750","Cre05.g235750","Cre06.g269950","Cre06.g250850","Cre06.g257800","Cre06.g285650","Cre06.g294200","Cre06.g295700","Cre07.g316850","Cre07.g325716","Cre07.g336650","Cre07.g347600","Cre07.g355200","Cre07.g312350","Cre07.g314650","Cre07.g338000","Cre07.g350550","Cre08.g366400","Cre08.g374050","Cre08.g368050","Cre09.g397327","Cre09.g397845")
length(dna.example)


test.gene <- "Cre01.g000150"
test.gene <- "Cre01.g000250"
subset(all.go, GID == test.gene)
subset(go.data.frame, GID == test.gene)
subset(symbol.data.frame, GID == test.gene)

library(GO.db)
"GO:0055085" %in% Lkeys(GO.db::GOTERM)
"GO:0046873" %in% Lkeys(GO.db::GOTERM)
"GO:0030001" %in% Lkeys(GO.db::GOTERM)
"GO:0016020" %in% Lkeys(GO.db::GOTERM)

head(all.go)


## Descargo de KEGG las proteinas anotadas con CHLREDRAFT 
## https://www.genome.jp/kegg-bin/show_organism?menu_type=genome_info&org=cre

## Generate TxDb package from gff3 file
library("GenomicFeatures")

## Generate chromosome.info data frame
library(seqinr)

cre.genome.data <- read.fasta(file = "Creinhardtii_281_v5.0.fa",seqtype = "DNA")
chromosome.names <- getName(cre.genome.data)
chromosome.lengths <- sapply(X=getSequence(cre.genome.data),FUN = length)

chromosomes.info <- data.frame(chrom=chromosome.names,length=chromosome.lengths,is_circular=FALSE)

## Meta data info
meta.data.info <- data.frame(name=c("Resource URL","Genome"),value=c("https://phytozome.jgi.doe.gov/","v5.0"))

cre.txdb <- makeTxDbFromGFF(file = "Creinhardtii_281_v5.6.gene_exons.gff3",format = "gff3",dataSource = "Phytozome",organism = "Chlamydomonas reinhardtii",taxonomyId = 3055,chrominfo = chromosomes.info,metadata = meta.data.info)

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

identity <- read.table(file="result_identity.txt",header = F,as.is = T)

get.first <- function(elto)
{
  return(elto[[1]])
}
cre.names <- sapply(strsplit(identity[[1]],split=".t"),get.first)

chlredraft.cre <- chlredraft[identity[[2]]]
names(chlredraft.cre) <- cre.names