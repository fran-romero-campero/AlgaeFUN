## Generate TxDb package from gtf file
library("GenomicFeatures")

## Generate chromosome.info data frame
library(seqinr)

genome.data <- read.fasta(file = "../genome/chromochloris_zofingiensis.fa",seqtype = "DNA")
chromosome.names <- getName(genome.data)
chromosome.lengths <- sapply(X=getSequence(genome.data),FUN = length)

chromosomes.info <- data.frame(chrom=chromosome.names,length=chromosome.lengths,is_circular=FALSE)

## Meta data info
meta.data.info <- data.frame(name=c("Resource URL","Genome"),
                             value=c("https://phytozome.jgi.doe.gov","v5.2.3.2"))

## Generate txdb object
microalgae.txdb <- makeTxDbFromGFF(file = "../annotation/chromochloris_zofingiensis.gtf",
                                   format = "gtf",
                                   dataSource = "Phytozome",
                                   organism = "Chromochloris zofingiensis",
                                   taxonomyId = 31302,
                                   chrominfo = chromosomes.info,
                                   metadata = meta.data.info)
head(as.data.frame(genes(microalgae.txdb)))

makeTxDbPackage(txdb = microalgae.txdb, 
                version = "0.1", 
                maintainer = "Francisco J. Romero-Campero <fran@us.es>", 
                author = "Pedro de los Reyes")

install.packages("./TxDb.Czofingiensis.Phytozome/", repos=NULL)
## loading packages
library(TxDb.Czofingiensis.Phytozome)
txdb <- TxDb.Czofingiensis.Phytozome

head(as.data.frame(genes(txdb)))
