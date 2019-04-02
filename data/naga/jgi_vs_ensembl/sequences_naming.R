library(seqinr)

naga.ensembl <- read.fasta(file="Nannochloropsis_gaditana_ensembl.fa",seqtype = "AA")
naga.ensembl.names.raw <- getAnnot(naga.ensembl)

naga.annot <- naga.ensembl.names.raw[[1]]

extract.ensembl.name <- function(naga.annot)
{
  return(strsplit(strsplit(x = naga.annot,split="gene:")[[1]][2],split=" ")[[1]][1])
}

naga.ensembl.names <- sapply(naga.ensembl.names.raw,FUN = extract.ensembl.name)

write.fasta(sequences = getSequence(naga.ensembl),names = naga.ensembl.names,file.out = "naga_ensembl.fa")

naga.jgi <- read.fasta(file = "Nannochloropsis_gaditana_jgi.fasta", seqtype = "AA")
naga.jgi.names.raw <- getAnnot(naga.jgi)

naga.jgi.names.raw[[1]]