## Author: Francisco J. Romero-Campero, Ana Belén Romero-Losada
## Date: November 2018
## Contact: fran@us.es

## Loading packages
library(ChIPseeker)
library(TxDb.Creinhardtii.Phytozome)
library(clusterProfiler)

## Auxiiary functions
extract.annotation <- function(annotation.str)
{
  return(strsplit(x = annotation.str,split=" \\(")[[1]][1])
}

input <- list(analysis = "genomic_locations", peaks="output_peaks.narrowPeak", promoter_length=1000,selected_genomic_features=c("Promoter","Gene Body"))
# gene_set

## Define transcript data base for the corresponding microalgae
txdb <- TxDb.Creinhardtii.Phytozome

## Load peak file
peaks <- readPeakFile(peakfile = input$peaks,header=FALSE)

## Define promoter region around TSS
promoter <- getPromoters(TxDb=txdb, upstream=input$promoter_length, downstream=input$promoter_length)

## Plot signal in promoter
plotAvgProf(tagMatrix, xlim=c(-input$promoter_length, input$promoter_length),
            xlab="Genomic Region (5'->3')", ylab = "Signal")


## Annotate genomic loci
peakAnno <- annotatePeak(peak, tssRegion=c(-input$promoter_length, input$promoter_length),
                         TxDb=txdb)#, annoDb="org.Creinhardtii.eg.db")

## Plot pie chart with annotation
plotAnnoPie(peakAnno)
#plotAnnoBar(peakAnno)
#vennpie(peakAnno)
#upsetplot(peakAnno)
#upsetplot(peakAnno, vennpie=TRUE)

## Plot distance to tss
plotDistToTSS(peakAnno,
              title="Distribution of genomic loci relative to TSS",ylab = "Genomic Loci (%) (5' -> 3')")

## Extract genes 
peak.annotation <- as.data.frame(peakAnno)
simple.annotation <- sapply(X = as.vector(peak.annotation$annotation),FUN = extract.annotation)
names(simple.annotation) <- NULL



genes.promoter <- peak.annotation$geneId[simple.annotation == "Promoter"]
genes.5utr <- peak.annotation$geneId[simple.annotation == "5' UTR"]
genes.3utr <- peak.annotation$geneId[simple.annotation == "3' UTR"]
genes.exon <- peak.annotation$geneId[simple.annotation == "Exon"]
genes.intron <- peak.annotation$geneId[simple.annotation == "Intron"]

selected.genes <- c(genes.promoter, genes.5utr,genes.3utr,genes.exon,genes.intron)
length(selected.genes)


unique(simple.annotation)



table(simple.annotation)/sum(table(simple.annotation))

subset(peak.annotation,annotation == "Promoter")$geneId





nrow(peak.annotation)

table(peak.annotation$annotation)


11316/nrow(peak.annotation)

promoter <- getPromoters(TxDb=txdb, upstream=input$promoter_upstream, downstream=input$promoter_downstream)
tagMatrixList <- lapply(input$peaks, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-input$promoter_upstream, input$promoter_downstream))

plotAvgProf(tagMatrixList, xlim=c(-1000, 1000), conf=0.95,resample=500, facet="row")

tagHeatmap(tagMatrixList, xlim=c(-1000, 1000), color=NULL)