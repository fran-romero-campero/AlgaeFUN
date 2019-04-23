## Generate TxDb package from gtf file
library("GenomicFeatures")

## Generate chromosome.info data frame
library(seqinr)

genome.data <- read.fasta(file = "../genome/coccomyxa_subellipsoidea.fa",seqtype = "DNA")
chromosome.names <- getName(genome.data)
chromosome.lengths <- sapply(X=getSequence(genome.data),FUN = length)

chromosomes.info <- data.frame(chrom=chromosome.names,length=chromosome.lengths,is_circular=FALSE)

## Meta data info
meta.data.info <- data.frame(name=c("Resource URL","Genome"),
                             value=c("https://phytozome.jgi.doe.gov/","v2.0"))

## Generate txdb object
microalgae.txdb <- makeTxDbFromGFF(file = "../annotation/coccomyxa_subellipsoidea.gtf",
                                   format = "gtf",
                                   dataSource = "Phytozome",
                                   organism = "Coccomyxa subellipsoidea",
                                   taxonomyId = 574566,
                                   chrominfo = chromosomes.info,
                                   metadata = meta.data.info)
genes(microalgae.txdb)

makeTxDbPackage(txdb = microalgae.txdb, 
                version = "0.1", 
                maintainer = "Francisco J. Romero-Campero <fran@us.es>", 
                author = "Ana B. Romero-Losada")

install.packages("./TxDb.Csubellipsoidea.Phytozome/", repos=NULL)
## loading packages
library(TxDb.Csubellipsoidea.Phytozome)
txdb <- TxDb.Csubellipsoidea.Phytozome

genes(txdb)
