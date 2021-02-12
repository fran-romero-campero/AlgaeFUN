
gene.link.function <- phytozome.gene.link
txdb <- TxDb.Creinhardtii.Phytozome
org.db <- org.Creinhardtii.eg.db

## Output table with gene annotation
annotations <- intersect(c("GO", "KO", "KOG", "ENZYME", "PANTHER","PFAM"),columns(org.db))
microalgae.annotation <- select(org.db,columns = annotations,keys=keys(org.db,keytype = "GID"))
genes.annotation <- microalgae.annotation

head(genes.annotation)
genes <- unique(genes.annotation$GID)
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
head(genes.annotation.links)

write.table(x = genes.annotation.download,file = "hlacustris_gene_annotation.tsv",sep="\t",row.names = FALSE,quote = F)
res <- read.table(file = "hlacustris_gene_annotation.tsv",sep="\t",header = T,as.is=T)
head(res)

write.table(x = genes.annotation.links,file = "hlacustris_gene_annotation_links.tsv",sep="\t",quote=FALSE)
res2 <- read.table(file = "hlacustris_gene_annotation_links.tsv",sep="\t",header = T,as.is=T)
head(res2)
