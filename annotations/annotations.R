library(TxDb.MpusillaCCMP1545.Phytozome)
library(org.MpusillaCCMP1545.eg.db)

split.commas <- function(annotation.str)
{
  return(strsplit(annotation.str,split=",")[[1]])
}

go.link <- function(go.term)
{
  link <- paste0("http://amigo.geneontology.org/amigo/term/", go.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           go.term, "</a>"),
                         collapse = "")
  return(complete.link)
}

## KO link
#https://www.genome.jp/dbget-bin/www_bget?ko:K00276
ko.link <- function(ko.term)
{
  link <- paste0("https://www.genome.jp/dbget-bin/www_bget?ko:", ko.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           ko.term, "</a>"),
                         collapse = "")
  return(complete.link)
}

## KOG link
## https://www.ncbi.nlm.nih.gov/Structure/cdd/KOG3720
kog.link <- function(kog.term)
{
  link <- paste0("https://www.ncbi.nlm.nih.gov/Structure/cdd/KOG3720", kog.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           kog.term, "</a>"),
                         collapse = "")
  return(complete.link)
}

## ENZYME link
## https://www.brenda-enzymes.org/enzyme.php?ecno=1.4.3.21
enzyme.link <- function(ec.term)
{
  link <- paste0("https://www.brenda-enzymes.org/enzyme.php?ecno=", ec.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           ec.term, "</a>"),
                         collapse = "")
  return(complete.link)
}

## PANTHER link
## http://www.pantherdb.org/panther/family.do?clsAccession=PTHR10638
panther.link <- function(panther.term)
{
  link <- paste0("http://www.pantherdb.org/panther/family.do?clsAccession=", panther.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           panther.term, "</a>"),
                         collapse = "")
  return(complete.link)
}

## PFAM link
## http://pfam.xfam.org/family/PF02728
pfam.link <- function(pfam.term)
{
  link <- paste0("http://pfam.xfam.org/family/", pfam.term)
  complete.link <- paste(c("<a href=\"",
                           link,
                           "\" target=\"_blank\">",
                           pfam.term, "</a>"),
                         collapse = "")
  return(complete.link)
}

no.link <- function(gene.name)
{
  return(gene.name)
}

gene.link.function <- no.link
txdb <- TxDb.MpusillaCCMP1545.Phytozome
org.db <- org.MpusillaCCMP1545.eg.db

## Output table with gene annotation
annotations <- intersect(c("GO", "KO", "KOG", "ENZYME", "PANTHER","PFAM"),columns(org.db))
microalgae.annotation <- select(org.db,columns = annotations,keys=keys(org.db,keytype = "GID"))
genes.annotation <- microalgae.annotation

head(genes.annotation)
genes <- unique(genes.annotation$GID)
length(genes)
genes.annotation.download <- data.frame(matrix(nrow=length(genes),ncol=(length(annotations)+1)))
genes.annotation.links <- genes.annotation.download 

colnames(genes.annotation.download)[1] <- "Gene ID"
colnames(genes.annotation.links)[1] <- "Gene ID"

for(i in 1:length(annotations))
{
  print(i)
  current.annotation <- annotations[i]
  
  colnames(genes.annotation.download)[i+1] <- current.annotation
  colnames(genes.annotation.links)[i+1] <- current.annotation
  
  if(current.annotation == "GO")
  {
    annotation.link <- go.link
  } else if (current.annotation == "KO")
  {
    annotation.link <- ko.link
  } else if (current.annotation == "KOG")
  {
    annotation.link <- kog.link
  } else if (current.annotation == "ENZYME")
  {
    annotation.link <- enzyme.link
  } else if (current.annotation == "PANTHER")
  {
    annotation.link <- panther.link
  } else if (current.annotation == "PFAM")
  {
    annotation.link <- pfam.link
  }
  
  for(j in 1:length(genes))
  {
    print(j)
    current.gene <- genes[j]
    genes.annotation.download[j,1] <- current.gene
    genes.annotation.links[j,1] <- gene.link.function(current.gene)
    current.gene.annotation <- sapply(unique(subset(genes.annotation, GID == current.gene)[[current.annotation]]),split.commas)
    if(is.na(current.gene.annotation[1]))
    {
      genes.annotation.download[j,(i+1)] <- ""
      genes.annotation.links[j,(i+1)] <- ""
    } else
    {
      genes.annotation.download[j,(i+1)] <- paste(current.gene.annotation,collapse=" ")
      genes.annotation.links[j,(i+1)] <- paste(sapply(current.gene.annotation,annotation.link),collapse=" ")
    }
  }
}

genes.annotation.download[genes.annotation.download == "character(0)"] <- ""
head(genes.annotation.download)

rownames(genes.annotation.links) <- genes.annotation.download$`Gene ID`
#head(genes.annotation.links)

write.table(x = genes.annotation.download,file = "mpusilla_gene_annotation.tsv",sep="\t",row.names = FALSE,quote = F)
res <- read.table(file = "mpusilla_gene_annotation.tsv",sep="\t",header = T,as.is=T,comment.char = "")
head(res)

write.table(x = genes.annotation.links,file = "mpusilla_gene_annotation_links.tsv",sep="\t",quote=FALSE)
res2 <- read.table(file = "mpusilla_gene_annotation_links.tsv",sep="\t",header = T,as.is=T,comment.char = "")
head(res2)
