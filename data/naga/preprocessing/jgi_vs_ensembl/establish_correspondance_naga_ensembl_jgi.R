## Authors: Francisco J. Romero-Campero
## Email: fran@us.es
## Date: March 2019

arguments <- commandArgs(trailingOnly = TRUE)
group <- as.numeric(arguments[1])

## Loading require packages seqinr and Biostrings
library(seqinr)
library(Biostrings)

## Loading substitution matrix PAM250
data(BLOSUM62)

## Loading details for source protein sequence
source.data <- read.fasta(file="naga_ensembl.fa", seqtype="AA")
source.names <- getName(source.data)
source.seqs <- getSequence(source.data)
n.source <- length(source.seqs)

## Loading details for target protein sequences  
target.data <- read.fasta(file="naga_jgi.fa", seqtype="AA")
target.names <- getName(target.data)
target.seqs <- getSequence(target.data)
n.target <- length(target.names)

## Initializing two vectors that will store the protein names of the 20 most similar 
## target sequences and their identity

size.batch <- 1000
pot.orthologs.identities <- vector(mode = "numeric",length = size.batch)
pot.orthologs.names <- vector(mode = "character",length = size.batch)

## Using a loop to go through all the target proteins    

lower.bound <- (group - 1)*size.batch + 1
upper.bound <- group * size.batch
k <- 1
for (i in lower.bound:upper.bound)
{
  print(i)
  source.protein.name <- source.names[[i]]
  source.protein.seq <- c2s(source.seqs[[i]])
  
  pot.orthologs.identity <- 0
  pot.orthologs.name <- NA
  
  j <- 1
  while(j <= n.target & pot.orthologs.identity != 1) #for(j in 1:n.target)  
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
    j <- j + 1
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
