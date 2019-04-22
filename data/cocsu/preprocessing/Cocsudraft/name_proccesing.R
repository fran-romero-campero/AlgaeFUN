library(seqinr)


seq.raw.data <- read.fasta(file = "GCF_000258705.1_Coccomyxa_subellipsoidae_v2.0_protein.faa", seqtype = "AA")
seqs <- getSequence(seq.raw.data)
seq.annot <- getAnnot(seq.raw.data)


  #extract.cocsu.keggid <- function (seq.annot)
 #{
    res<- c()
  i<-1
    for (i in 1:length(seq.annot))
  {
      current.annot <- seq.annot[[i]]
      
    if (grepl(pattern = "COCSU", x = current.annot))
    {
      if (grepl(pattern = "partial", x = current.annot))
      {
      splited.annot <- strsplit(x = current.annot, split = "COCSU")[[1]][2]
      splited.annot.2 <- strsplit(x = splited.annot, split = " ")[[1]][1]
      splited.annot.3 <- strsplit(x = splited.annot.2, split= "," )[[1]][1]
      res[i] <- paste("COCSU", splited.annot.3, sep = "")
      
      } else
      splited.annot.j <- strsplit(x = current.annot, split = "COCSU")[[1]][2]
      splited.annot.g <- strsplit(x = splited.annot.j, split = " ")[[1]][1]
      res[i] <- paste("COCSU", splited.annot.g, sep = "")
    } else 
    {
    res[i] <- NA
    }
  
    #return(res)
  }
#}

res[34:999]
length(res)
length(seq.annot)

write.fasta(sequences = seqs, names = res, file.out = "cocsudraft.fa")
