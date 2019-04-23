## Generate TxDb package from gtf file
library("GenomicFeatures")

## Generate chromosome.info data frame
library(seqinr)

genome.data <- read.fasta(file = "../genome/phaeodactylum_tricornutum.fa",seqtype = "DNA")
chromosome.names <- getName(genome.data)
chromosome.lengths <- sapply(X=getSequence(genome.data),FUN = length)

chromosomes.info <- data.frame(chrom=chromosome.names,length=chromosome.lengths,is_circular=FALSE)

## Meta data info
meta.data.info <- data.frame(name=c("Resource URL","Genome"),
                             value=c("https://protists.ensembl.org/","v2.0"))

## Generate txdb object
microalgae.txdb <- makeTxDbFromGFF(file = "../annotation/phaeodactylum_tricornutum.gtf",
                                   format = "gtf",
                                   dataSource = "Ensembl Protists",
                                   organism = "Phaeodactylum tricornutum",
                                   taxonomyId = 556484,
                                   chrominfo = chromosomes.info,
                                   metadata = meta.data.info)


makeTxDbPackage(txdb = microalgae.txdb, 
                version = "0.1", 
                maintainer = "Francisco J. Romero-Campero <fran@us.es>", 
                author = "Francisco J. Romero-Campero")

install.packages("./TxDb.Ptricornutum.Ensembl.Protists/", repos=NULL)
## loading packages
library(TxDb.Ptricornutum.Ensembl.Protists)
txdb <- TxDb.Ptricornutum.Ensembl.Protists

txdb
