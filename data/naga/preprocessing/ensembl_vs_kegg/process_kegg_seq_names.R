library(seqinr)

seq.raw.data <- read.fasta(file = "naga_kegg_raw.fa",seqtype = "AA")
seq.names <- getName(seq.raw.data)
seqs <- getSequence(seq.raw.data)

write.fasta(sequences = seqs,names = seq.names,file.out = "naga_kegg.fa")
