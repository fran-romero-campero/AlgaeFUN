## Authors: Francisco J. Romero-Campero
## Email: fran@us.es
## Date: December 2018

arguments <- commandArgs(trailingOnly = TRUE)
group <- as.numeric(arguments[1])

## Loading require packages seqinr and Biostrings
library(seqinr)
library(Biostrings)

## Loading substitution matrix PAM250
data(BLOSUM62)

## Loading details for source protein sequence
source.data <- read.fasta(file="Creinhardtii_281_v5.6.protein.fa", seqtype="AA")
source.names <- getName(source.data)
source.seqs <- getSequence(source.data)
n.source <- length(source.seqs)

## Loading details for target protein sequences  
target.data <- read.fasta(file="GCF_000002595.1_v3.0_translated_cds.faa", seqtype="AA")
target.names <- getName(target.data)
target.seqs <- getSequence(target.data)
n.target <- length(target.names)

## Initializing two vectors that will store the protein names of the 20 most similar 
## target sequences and their identity

pot.orthologs.identities <- vector(mode = "numeric",length = n.source)
pot.orthologs.names <- vector(mode = "character",length = n.source)

## Using a loop to go through all the target proteins    

lower.bound <- (group - 1)*1000 + 1
upper.bound <- group * 1000
k <- 1
for (i in lower.bound:upper.bound)
{
  print(i)
  source.protein.name <- source.names[[i]]
  source.protein.seq <- c2s(source.seqs[[i]])
  
  pot.orthologs.identity <- 0
  pot.orthologs.name <- NA
  
  for(j in 1:n.target)  
  {
    target.protein.name <- target.names[[j]]
    target.protein.seq <- c2s(target.seqs[[j]])
    
    ## Computing best global aligment using Needleman-Wunsch algorithm and the identity between the 
    ## current target protein and the source protein
    res.align <- pairwiseAlignment(pattern=source.protein.seq, subject=target.protein.seq, substitutionMatrix=BLOSUM62, gapOpening=0)
    source.protein.align <- s2c(as.character(pattern(res.align)))
    target.protein.align <- s2c(as.character(subject(res.align)))
    identity <- sum(target.protein.align == source.protein.align)/length(source.protein.align)
    
    if(identity > pot.orthologs.identity)
    {
      pot.orthologs.identity <- identity
      pot.orthologs.name <- target.protein.name
    }
  }
  
  pot.orthologs.identities[k] <- pot.orthologs.identity
  pot.orthologs.names[k] <- pot.orthologs.name
  print(pot.orthologs.identity)
  print("----------")
  k <- k + 1
}

res.identity <- data.frame(source.names[lower.bound:upper.bound],
                           pot.orthologs.names,
                           pot.orthologs.identities)

write.table(x = res.identity,file = paste(c("result_identity_",group,".txt"),collapse=""),row.names=FALSE,col.names=FALSE)



setdiff(target.seqs[[j]],colnames(BLOSUM62))
colnames(BLOSUM62)
