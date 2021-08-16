## Generate example of genomic loci
library(TxDb.Smuscicola.pub)
txdb <- TxDb.Smuscicola.pub

nrow(as.data.frame(cds(txdb)))

number.random.peaks <- 500

genes.coordinates <- as.data.frame(genes(txdb))
unique(genes.coordinates$seqnames)
nrow(genes.coordinates)
genes.coordinates <- subset(genes.coordinates, !(seqnames %in% c("Pt","Mt")))
nrow(genes.coordinates)

selected.id <- sample(x = 1:nrow(genes.coordinates),size = number.random.peaks,replace = F)

selected.genes.coordinates <- genes.coordinates[selected.id,]
head(selected.genes.coordinates)

distance.tss <- sample(x = 1:20,size = number.random.peaks,replace = TRUE)*100
up.down.stream <- sample(x = c(-1,1),size = number.random.peaks,replace = TRUE)
region.width <- sample(x = 1:3,size = number.random.peaks,replace = TRUE)*100

selected.genes.coordinates$start <- selected.genes.coordinates$start + up.down.stream * distance.tss

random.peaks <- data.frame(selected.genes.coordinates$seqnames, 
                           selected.genes.coordinates$start - region.width/2,
                           selected.genes.coordinates$start + region.width/2)

colnames(random.peaks) <- NULL
head(random.peaks)

random.peaks <- random.peaks[(random.peaks[,2] > 0 & random.peaks[,3] > 0),]

write.table(x = random.peaks,file = "unsorted_example_genomic_regions_smuscicola.txt",quote = F,sep = "\t",row.names = F, col.names = F)



genomic.regions <- readPeakFile(peakfile = "example_files/example_genomic_regions_mpusilla.txt")
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, 
                         downstream=1000)

peakAnno <- annotatePeak(peak = genomic.regions, 
                         tssRegion=c(-1000, 1000),
                         TxDb=txdb)

genes.txdb <- genes(txdb)
head(genes.txdb)
genes.names <- names(genes.txdb)


genes.data <- as.data.frame(genes(txdb))
genes.1 <- genes.data$gene_id

genes.data.2 <- as.data.frame(genes(txdb,single.strand.genes.only=FALSE))
head(genes.data.2)
genes.2 <- genes.data.2$group_name

setdiff(genes.2,genes.1)
