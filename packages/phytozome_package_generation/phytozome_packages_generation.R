## Author: Francisco J. Romero-Campero, Ana Belén Romero-Losada
## Date: November 2018
## Contact: fran@us.es

## Documentation on how to generate org.db packages can be found on these links:
## http://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html
## 

## Replace ' with prime in the info.txt file

microalgae.info <- read.table(file="MpusillaCCMP1545_228_v3.0.annotation_info.txt",header=T,comment.char = "",sep = "\t",as.is=T)
head(microalgae.info)

write.table(x = unique(microalgae.info$locusName),file = "mpusilla_universe.txt",quote = F,row.names = F,col.names = F)
length(unique(microalgae.info$locusName))

## Generate and write output data frame
id.microalgae.name <- data.frame(GID=microalgae.info$locusName,GENENAME=microalgae.info$locusName,stringsAsFactors = FALSE)
head(id.microalgae.name)

## Generate SYMBOL data frame
symbol.data.frame <- data.frame(GID=microalgae.info$locusName,SYMBOL=microalgae.info$locusName,stringsAsFactors = FALSE)
head(symbol.data.frame)
symbol.data.frame[1:20,]

## Remove duplicated rows
nrow(symbol.data.frame)
symbol.data.frame <- symbol.data.frame[!duplicated(symbol.data.frame),]
nrow(symbol.data.frame)


## Generate GO data.frame to create org.Db package
## Before reading the GO annotation file I have to replace the symbol ' by the word prime
## to avoid the error "EOF within quoted string"

## According to the vignettes of AnnotationForge bioconductor package for "Making Organisms
## Packages: "However to use the goTable argument, you have to follow a strict convention with 
## the data. Such a data.frame must have three columns only and these must correspond to the 
## gene id, GO id and evidence codes. These columns also have to be named as “GID”, “GO” and “EVIDENCE”

gid <- microalgae.info$locusName
go <- microalgae.info$GO
evidence <- rep("ISS",length(gid))

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
nrow(go.data.frame)
go.data.frame <- go.data.frame[!duplicated(go.data.frame),]
go.data.frame[1:40,]
nrow(go.data.frame)
length(unique(go.data.frame$GID))

## For the rest of the data frames. According to the help info "the 1st column of every 
## data.frame must be labeled GID, and correspond to a gene ID that is universal for the 
## entire set of data.frames.  The GID is how the different tables will be joined internally"

## Generate enzyme data.frame

enzyme <- microalgae.info$KEGG.ec
 
enzyme.data.frame <- data.frame(GID=gid,ENZYME=enzyme,stringsAsFactors = FALSE)
head(enzyme.data.frame)

## Remove duplicated rows
nrow(enzyme.data.frame)
enzyme.data.frame <- enzyme.data.frame[!duplicated(enzyme.data.frame),]
head(enzyme.data.frame)
nrow(enzyme.data.frame)
length(unique(enzyme.data.frame$GID))

## Generate KOG data.frame
kog <- microalgae.info$KOG

kog.data.frame <- data.frame(GID=gid,KOG=kog,stringsAsFactors = FALSE)
head(kog.data.frame)

## Remove duplicated rows
kog.data.frame <- kog.data.frame[!duplicated(kog.data.frame),]

## Generate KO data.frame
ko <- microalgae.info$KO

ko.data.frame <- data.frame(GID=gid,KO=ko,stringsAsFactors = FALSE)
head(ko.data.frame)

## Remove duplicated rows
ko.data.frame <- ko.data.frame[!duplicated(ko.data.frame),]
head(ko.data.frame)


## Generate PFAM data.frame
pfam <- microalgae.info$Pfam

pfam.data.frame <- data.frame(GID=gid,PFAM=pfam,stringsAsFactors = FALSE)
head(pfam.data.frame)

## Remove duplicated rows
pfam.data.frame <- pfam.data.frame[!duplicated(pfam.data.frame),]
head(pfam.data.frame)

## Generate PANTHER data.frame
panther <- microalgae.info$Panther

panther.data.frame <- data.frame(GID=gid,PANTHER=panther,stringsAsFactors = FALSE)
head(panther.data.frame)

## Remove duplicated rows
panther.data.frame <- panther.data.frame[!duplicated(panther.data.frame),]
head(panther.data.frame)

## Microalgae Taxonomy ID: 564608
## Go to https://www.ncbi.nlm.nih.gov/ in the search box type the name of your
## microalgae. Scroll down and go to Taxonomy at the very bottom. Click on Taxonomy.
## Click on your microalgae name and copy Taxonomy ID:

## Load require package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationForge")
library(AnnotationForge)

makeOrgPackage(go=go.data.frame,
               SYMBOL=symbol.data.frame,
               ENZYME=enzyme.data.frame,
               KOG=kog.data.frame,
               KO=ko.data.frame,
               PFAM=pfam.data.frame,
               PANTHER=panther.data.frame,
               version = "0.1",
               maintainer = "Francisco J. Romero-Campero <fran@us.es>",
               author = "Christina Arvanitidou",
               outputDir = ".", 
               tax_id = "564608",
               genus = "Micromonas",
               species = "pusillaCCMP1545",
               goTable = "go",
               verbose = TRUE)

install.packages("./org.MpusillaCCMP1545.eg.db/", repos=NULL)
remove.packages("org.Creinhardtii.eg.db")


library(org.MpusillaCCMP1545.eg.db)
columns(org.MpusillaCCMP1545.eg.db)
head(select(org.MpusillaCCMP1545.eg.db,columns = c("GO"),keys=keys(org.MpusillaCCMP1545.eg.db,keytype = "GID")))
head(select(org.MpusillaCCMP1545.eg.db,columns = c("GOALL"),keys=keys(org.MpusillaCCMP1545.eg.db,keytype = "GID")))
head(select(org.MpusillaCCMP1545.eg.db,columns = c("ENZYME"),keys=keys(org.MpusillaCCMP1545.eg.db,keytype = "GID")))
head(select(org.MpusillaCCMP1545.eg.db,columns = c("EVIDENCE"),keys=keys(org.MpusillaCCMP1545.eg.db,keytype = "GID")))
head(select(org.MpusillaCCMP1545.eg.db,columns = c("EVIDENCEALL"),keys=keys(org.MpusillaCCMP1545.eg.db,keytype = "GID")))
head(select(org.MpusillaCCMP1545.eg.db,columns = c("KOG"),keys=keys(org.MpusillaCCMP1545.eg.db,keytype = "GID")))
head(select(org.MpusillaCCMP1545.eg.db,columns = c("ONTOLOGY"),keys=keys(org.MpusillaCCMP1545.eg.db,keytype = "GID")))
head(select(org.MpusillaCCMP1545.eg.db,columns = c("ONTOLOGYALL"),keys=keys(org.MpusillaCCMP1545.eg.db,keytype = "GID")))
head(select(org.MpusillaCCMP1545.eg.db,columns = c("SYMBOL"),keys=keys(org.MpusillaCCMP1545.eg.db,keytype = "GID")))

## Preprocess gff3 to generate gtf with gene_id and transcript_id 

## Corresponde gene id pacid
pacId <- microalgae.info$pacId
geneId <- microalgae.info$locusName
names(geneId) <- pacId


microalgae.gff3 <- read.table(file="MpusillaCCMP1545_228_v3.0.gene_exons.gff3",header=F,quote = "#",as.is=T)
microalgae.gtf <- microalgae.gff3
head(microalgae.gff3)

unique(microalgae.gff3$V3)

for(i in 1:nrow(microalgae.gff3))
{
  current.attributes <- strsplit(microalgae.gff3$V9[i],split=";")[[1]]
  if(microalgae.gff3$V3[i] == "gene")
  {
    gene.id <- strsplit(current.attributes[2],split="=")[[1]][2]
    microalgae.gtf$V9[i] <- paste("gene_id", paste("\"",gene.id,"\";",sep=""))
  } else if(microalgae.gff3$V3[i] == "mRNA")
  {
    gene.name <- strsplit(current.attributes[5],split="=")[[1]][2]
    gene.id <- substr(gene.name,start = 1,stop = nchar(gene.name)-8)
    transcript.id <- strsplit(current.attributes[2],"=")[[1]][2]
    microalgae.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep="")), collapse = " ")
  } else if(microalgae.gff3$V3[i] == "exon")
  {
    gene.id <- geneId[strsplit(current.attributes[3],split="=")[[1]][2]]
    transcript.name <- strsplit(current.attributes[2],split="=")[[1]][2]
    trancript.id <- substr(transcript.name,start = 1,stop = nchar(transcript.name) - 8)
    exon.number <- strsplit(strsplit(current.attributes[1],split="=")[[1]][2],split="exon.")[[1]][2]
    microalgae.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep=""),"exon_number",paste("\"",exon.number,"\";",sep="")), collapse = " ")
  } else if(microalgae.gff3$V3[i] == "five_prime_UTR" || microalgae.gff3$V3[i] == "three_prime_UTR" || microalgae.gff3$V3[i] == "CDS")
  {
    gene.id <- geneId[strsplit(current.attributes[3],split="=")[[1]][2]]
    transcript.name <- strsplit(current.attributes[2],split="=")[[1]][2]
    trancript.id <- substr(transcript.name,start = 1,stop = nchar(transcript.name) - 8)
    microalgae.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep="")), collapse = " ")
  }
}

head(microalgae.gtf)

write.table(x = microalgae.gtf,file = "micromonas_pusillaCCMP1545.gtf",sep = "\t",row.names = F,col.names = F,quote = F)

## Generate TxDb package from gff3 file
library("GenomicFeatures")

## Generate chromosome.info data frame
library(seqinr)

microalgae.genome.data <- read.fasta(file = "../../genomes/mpusillaCCMP1545.fa",seqtype = "DNA")
chromosome.names <- getName(microalgae.genome.data)
chromosome.lengths <- sapply(X=getSequence(microalgae.genome.data),FUN = length)

chromosomes.info <- data.frame(chrom=chromosome.names,length=chromosome.lengths,is_circular=FALSE)

## Meta data info
meta.data.info <- data.frame(name=c("Resource URL","Genome"),value=c("https://phytozome.jgi.doe.gov/","v3.0"))

microalgae.txdb <- makeTxDbFromGFF(file = "micromonas_pusillaCCMP1545.gtf",format = "gtf",
                                   dataSource = "Phytozome",organism = "Micromonas pusillaCCMP1545",
                                   taxonomyId = 564608,chrominfo = chromosomes.info,metadata = meta.data.info)

microalgae.txdb
genes(microalgae.txdb)

?makeTxDbPackage

makeTxDbPackage(txdb = microalgae.txdb, version = "0.1", maintainer = "Francisco J. Romero-Campero <fran@us.es>", 
                author = "Christina Arvanitidou")

install.packages("./TxDb.MpusillaCCMP1545.Phytozome/", repos=NULL)