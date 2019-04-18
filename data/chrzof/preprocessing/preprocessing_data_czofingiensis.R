## Before reading the info file # was removed and ' were replaced by prime
czof.data <- read.table(file = "Czofingiensis_461_v5.2.3.2.annotation_info.txt",header = T,sep = "\t",as.is = T,fill = T)

## Check that the data frame has the same number of rows as the corresponding file
nrow(czof.data)

## Generate GENENAME data frame
gene.name.data.frame <- data.frame(GID=czof.data$locusName,GENENAME=czof.data$locusName,stringsAsFactors = FALSE)
nrow(gene.name.data.frame)
head(id.czof.name)
gene.name.data.frame <- gene.name.data.frame[!duplicated(gene.name.data.frame),]
nrow(gene.name.data.frame)

## Generate SYMBOL data frame
symbol.data.frame <- data.frame(GID=czof.data$locusName,GENENAME=czof.data$locusName,stringsAsFactors = FALSE)
head(symbol.data.frame)
nrow(symbol.data.frame)
symbol.data.frame <- symbol.data.frame[!duplicated(symbol.data.frame),]
nrow(symbol.data.frame)

## Generate PFAM data frame
gid <- czof.data$locusName
pfam <- czof.data$Pfam

pfam.data.frame <- data.frame(GID=gid,PFAM=pfam,stringsAsFactors = FALSE)
head(pfam.data.frame)

pfam.data.frame <- subset(pfam.data.frame, PFAM != "")
head(pfam.data.frame)

pfam.expanded.data.frame <- c()
for(i in 1:nrow(pfam.data.frame))
{
  pfam.splitted <- strsplit(pfam.data.frame$PFAM[i],split=",")[[1]]
  
  if(length(pfam.splitted) > 0)
  {
    for(j in 1:length(pfam.splitted))
    {
      pfam.expanded.data.frame <- rbind(pfam.expanded.data.frame,c(pfam.data.frame[i,1],pfam.splitted[j]))
    }
  } else
  {
    pfam.expanded.data.frame <- rbind(pfam.expanded.data.frame,c(pfam.data.frame[i,1],""))
  }
}

head(pfam.expanded.data.frame)

pfam.data.frame <- data.frame(GID=pfam.expanded.data.frame[,1],
                              PFAM=pfam.expanded.data.frame[,2],
                              stringsAsFactors = FALSE)

head(pfam.data.frame)
pfam.data.frame <- pfam.data.frame[!duplicated(pfam.data.frame),]
pfam.data.frame[1:40,]
nrow(pfam.data.frame)
length(unique(pfam.data.frame$GID))

## Generate PANTHER data frame
gid <- czof.data$locusName
panther <- czof.data$Panther

panther.data.frame <- data.frame(GID=gid,PANTHER=panther,stringsAsFactors = FALSE)
head(panther.data.frame)

panther.data.frame <- subset(panther.data.frame, PANTHER != "")
head(panther.data.frame)

panther.expanded.data.frame <- c()
for(i in 1:nrow(panther.data.frame))
{
  panther.splitted <- strsplit(panther.data.frame$PANTHER[i],split=",")[[1]]
  
  if(length(panther.splitted) > 0)
  {
    for(j in 1:length(panther.splitted))
    {
      panther.expanded.data.frame <- rbind(panther.expanded.data.frame,
                                           c(panther.data.frame[i,1],panther.splitted[j]))
    }
  } else
  {
    panther.expanded.data.frame <- rbind(panther.expanded.data.frame,c(panther.data.frame[i,1],""))
  }
}

head(panther.expanded.data.frame)

panther.data.frame <- data.frame(GID=panther.expanded.data.frame[,1],
                                 PANTHER=panther.expanded.data.frame[,2],
                                 stringsAsFactors = FALSE)

head(panther.data.frame)
panther.data.frame <- panther.data.frame[!duplicated(panther.data.frame),]
panther.data.frame[1:40,]
nrow(panther.data.frame)
length(unique(panther.data.frame$GID))

## Generate KOG data frame
kog.data.frame <- data.frame(GID=czof.data$locusName,KOG=czof.data$KOG,stringsAsFactors = FALSE)
nrow(kog.data.frame)
head(kog.data.frame)
kog.data.frame <- subset(kog.data.frame, KOG != "")
head(kog.data.frame)
nrow(kog.data.frame)
kog.data.frame <- kog.data.frame[!duplicated(kog.data.frame),]
nrow(kog.data.frame)

## Generate EC data frame
ec.data.frame <- data.frame(GID=czof.data$locusName,EC=czof.data$ec,stringsAsFactors = FALSE)
nrow(ec.data.frame)
head(ec.data.frame)
ec.data.frame <- subset(ec.data.frame, EC != "")
head(ec.data.frame)
nrow(ec.data.frame)
ec.data.frame <- ec.data.frame[!duplicated(ec.data.frame),]
nrow(ec.data.frame)

## Generate KO data frame
ko.data.frame <- data.frame(GID=czof.data$locusName,KO=czof.data$KO,stringsAsFactors = FALSE)
nrow(ko.data.frame)
head(ko.data.frame)
ko.data.frame <- subset(ko.data.frame, KO != "")
head(ko.data.frame)
nrow(ko.data.frame)
ko.data.frame <- ko.data.frame[!duplicated(ko.data.frame),]
nrow(ko.data.frame)

## Generate GO data frame
gid <- czof.data$locusName
go <- czof.data$GO

go.data.frame <- data.frame(GID=gid,GO=go,stringsAsFactors = FALSE)
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
                            EVIDENCE=rep("ISS",nrow(go.expanded.data.frame)),
                            stringsAsFactors = FALSE)

head(go.data.frame)
go.data.frame <- go.data.frame[!duplicated(go.data.frame),]
go.data.frame[1:40,]
nrow(go.data.frame)
length(unique(go.data.frame$GID))

## Chromochloris zofigiensis Taxonomy ID: 31302

## Load require package
library(AnnotationForge)

makeOrgPackage(go=go.data.frame,
               SYMBOL=symbol.data.frame,
               ENZYME=ec.data.frame,
               KOG=kog.data.frame,
               KO=ko.data.frame,
               PANTHER=panther.data.frame,
               PFAM=pfam.data.frame,
               version = "0.1",
               maintainer = "Francisco J. Romero-Campero <fran@us.es>",
               author = "Christina Arvanitidou",
               outputDir = ".", 
               tax_id = "31302",
               genus = "Chromochloris",
               species = "zofingiensis",
               goTable = "go",
               verbose = TRUE)

install.packages("./org.Czofingiensis.eg.db/", repos=NULL)

library(org.Czofingiensis.eg.db)
columns(org.Czofingiensis.eg.db)
head(select(org.Czofingiensis.eg.db,columns = c("EC"),keys=keys(org.Czofingiensis.eg.db,keytype = "GID")))
head(select(org.Czofingiensis.eg.db,columns = c("GO"),keys=keys(org.Czofingiensis.eg.db,keytype = "GID")))
head(select(org.Czofingiensis.eg.db,columns = c("KOG"),keys=keys(org.Czofingiensis.eg.db,keytype = "GID")))
head(select(org.Czofingiensis.eg.db,columns = c("KO"),keys=keys(org.Czofingiensis.eg.db,keytype = "GID")))

## Preprocess gff3 to generate gtf with gene_id and transcript_id 
czof.gff3 <- read.table(file="Czofingiensis_461_v5.2.3.2.gene_exons.gff3",header=F,quote = "#",as.is=T)
czof.gtf <- czof.gff3
head(czof.gff3)

unique(czof.gff3$V3)

for(i in 1:nrow(czof.gff3))
{
  current.attributes <- strsplit(czof.gff3$V9[i],split=";")[[1]]
  if(czof.gff3$V3[i] == "gene")
  {
    gene.id <- strsplit(current.attributes[2],split="=")[[1]][2]
    czof.gtf$V9[i] <- paste("gene_id", paste("\"",gene.id,"\";",sep=""))
  } else if(czof.gff3$V3[i] == "mRNA")
  {
    transcript.id <- strsplit(current.attributes[2],"=")[[1]][2]
    gene.id <- strsplit(transcript.id,split="\\.")[[1]][1]
    czof.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep="")), collapse = " ")
  } else if(czof.gff3$V3[i] == "exon")
  {
    #gene.id <- substr(strsplit(current.attributes[1],split="=")[[1]][2],start = 1,stop = 10)
    trancript.id <- substr(strsplit(current.attributes[1],split="=")[[1]][2],start = 1,stop = 13)
    gene.id <- strsplit(transcript.id,split="\\.")[[1]][1]
    exon.number <- strsplit(strsplit(current.attributes[1],split="=")[[1]][2],split="exon.")[[1]][2]
    czof.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep=""),"exon_number",paste("\"",exon.number,"\";",sep="")), collapse = " ")
  } else if(czof.gff3$V3[i] == "CDS")
  {
    #gene.id <- substr(strsplit(current.attributes[2],split="=")[[1]][2],start = 1,stop = 10)
    transcript.id <- substr(strsplit(current.attributes[2],split="=")[[1]][2],start = 1,stop = 13)
    gene.id <- strsplit(transcript.id,split="\\.")[[1]][1]
    czof.gtf$V9[i] <- paste(c("gene_id",paste("\"",gene.id,"\";",sep=""),"transcript_id",paste("\"",transcript.id,"\";",sep="")), collapse = " ")
  }
}

head(czof.gtf)

write.table(x = czof.gtf,file = "chromochloris_zofingiensis.gtf",sep = "\t",row.names = F,col.names = F,quote = F)









