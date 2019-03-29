library(seqinr)

raw.data <- read.fasta(file = "phatri_draft_raw.fa",seqtype = "AA")
seq.names <- getAnnot(raw.data)
seq.names

seqs <- getSequence(raw.data)

i <- 1

gene.annot <- seq.names[[i]]

extract.draft <- function(gene.annot)
{
  return(strsplit(strsplit(gene.annot, split="locus_tag=")[[1]][2],split="]")[[1]][1])  
}

extract.draft(gene.annot)

new.names <- sapply(seq.names, extract.draft)

write.fasta(sequences = seqs,names = new.names,file.out = "phatri_proteins_draft.fa")

phatri.ensembl <- read.fasta(file = "Phaeodactylum_tricornutum_pep.fa",seqtype = "AA")
phatri.seqs <- getSequence(phatri.ensembl)
phatri.annot <- getAnnot(phatri.ensembl)

extract.phatri <- function(gene.annot)
{
  return(strsplit(strsplit(gene.annot,split = "gene:")[[1]][2],split=" ")[[1]][1])
}

phatri.names <- sapply(X = phatri.annot, FUN = extract.phatri)

write.fasta(sequences = phatri.seqs, names = phatri.names,file.out = "phatri_proteins_ensembl.fa")
