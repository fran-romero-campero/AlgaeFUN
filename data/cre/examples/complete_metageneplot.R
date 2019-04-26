## Author: Francisco J. Romero-Campero
## Email: fran@us.es

## This script receives as input files in bigWig format and generate metageneplots. 
## It assumes that the input files are located in the same folder as this file 
## (complete_metageneplot.R) so you can set the working directory to the source file
## location (In Rstudio menus: Session -> Set Working Directoy -> To Source File Location)

## This script uses the libraries ChIPpeakAnno, rtracklayer
## and TxDb.Athaliana.BioMart.plantsmart28. It is best executed in RStudio.
library(ChIPpeakAnno)
library(rtracklayer)
library(TxDb.Athaliana.BioMart.plantsmart28)

## If these packages are not installed try:
## source("https://bioconductor.org/biocLite.R")
## biocLite("ChIPpeakAnno")
## biocLite("rtracklayer")
## biocLite("TxDb.Athaliana.BioMart.plantsmart28")

## The genomic features of Arabidopsis thaliana obtained from the plants mart 28
## are used. 
txdb <- TxDb.Athaliana.BioMart.plantsmart28

## The different sets of genes are read. 
only.h2.repressed <- as.vector(read.table(file="only_h2_repressed.txt")[[1]])
length(only.h2.repressed)
only.h2.active <- as.vector(read.table(file="only_h2_active.txt")[[1]])
length(only.h2.active)
h2.k27 <- as.vector(read.table(file="h2_k27.txt")[[1]])
length(h2.k27)
only.k27 <- as.vector(read.table(file="only_k27.txt")[[1]])
length(only.k27)

## Specification of the set of genes to consider to construct the metageneplot
genes <- c(only.k27)
#genes <- c(h2.k27) ## Example for the set of genes with both marks H2AK121ub and H3K27me3
length(genes)

## Extraction of the genomic features of the specified genes.
genes.data <- subset(genes(txdb,columns=c("tx_id", "tx_name","gene_id")), gene_id %in% genes)
head(genes.data)

## Extraction of the TSS 
genes.tss <- resize(genes.data, width=1, fix='start')
head(genes.tss)

## Centering around TSS 2kb upstream and downstream
around.genes.tss <- genes.tss
start(around.genes.tss) <- start(genes.tss) - 2000
end(around.genes.tss) <- end(genes.tss) + 2000
head(around.genes.tss)
feature.recentered <- around.genes.tss

## BigWig files for the two conditions to compare
bigwig.files <- c("col_0_h3k27_rep1.bw",
                  "clfswn_h3k27_rep1.bw")

#bigwig.files <- c("col_0_h2aub_7dag_2.bw",     ## Example for comparing triple mutant 
#                  "atbmi1abc_h2aub_7dag_2.bw") ## AtBMI1abc with Col-0


## Importing bigWig files and naming data (you may want to change the 
## names).
cvglists <- sapply(bigwig.files, import, 
                   format="BigWig", 
                   which=feature.recentered, 
                   as="RleList")
names(cvglists) <- c("Col-0", "clfswn")

#names(cvglists) <- c("Col-0", "atbmi1abc") ## Example for comparing atbmi1abc with Col-0


## Extracting the signal 2kb from the center with 50 tile (This may be slow)
sig <- featureAlignedSignal(cvglists, around.genes.tss, 
                            upstream=2000, downstream=2000,n.tile=50) 

## Plotting the results around TSS (you may have to change the names of the conditions)
plot(colMeans(sig[["Col-0"]],na.rm = TRUE),type="l",col="blue",lwd=3,ylab="",cex.lab=2,axes=FALSE,xlab="",ylim=c(0,830))
lines(colMeans(sig[["clfswn"]],na.rm = TRUE),type="l",col="red",lwd=3)
axis(side = 1,labels = c("-2000","-1000","TSS","1000","2000"),at = c(1,13,26,39,50),lwd=2,cex=1.5,las=2,cex=2)
axis(side= 2, labels = seq(from=0,to=1000,by=200),at=seq(from=0,to=1000,by=200),lwd=2,cex=1.5)
mtext(text = "RPKM",side = 2,cex=2,line = 2.3,family="bold")
legend("topright",legend = c("Col-0",expression(italic("clf/swn"))),col = c("blue","red"),lwd=3,cex=1)

## example plot for atbmi1abc vs Col-0
#plot(colMeans(sig[["Col-0"]],na.rm = TRUE),type="l",col="blue",lwd=3,ylab="",cex.lab=2,axes=FALSE,xlab="",ylim=c(0,200))
#lines(colMeans(sig[["atbmi1abc"]],na.rm = TRUE),type="l",col="red",lwd=3)
#axis(side = 1,labels = c("-2000","-1000","TSS","1000","2000"),at = c(1,13,26,39,50),lwd=2,cex=1.5,las=2,cex=2)
#axis(side= 2, labels = seq(from=0,to=200,by=40),at=seq(from=0,to=200,by=40),lwd=2,cex=1.5)
#mtext(text = "RPKM",side = 2,cex=2,line = 2.3,family="bold")
#legend("topright",legend = c("Col-0",expression(italic("atbmi1abc"))),col = c("blue","red"),lwd=3,cex=1)

## Storing the signal around tss
signal.around.tss <- sig

## Extraction of the TES
genes.tes <- resize(genes.data, width=1, fix='end')
head(genes.tes)

## Centering around TES 2kb upstream and downstream
around.genes.tes <- genes.tes
start(around.genes.tes) <- start(genes.tes) - 2000
end(around.genes.tes) <- end(genes.tes) + 2000
head(around.genes.tes)

## Extracting the signal 2kb from the center with 50 tile (This may be slow)
sig <- featureAlignedSignal(cvglists, around.genes.tes, 
                            upstream=2000, downstream=2000,n.tile=50) 

## Plotting the results around TES (you may have to change the names of the conditions)
plot(colMeans(sig[["Col-0"]],na.rm = TRUE),type="l",col="blue",lwd=3,ylab="",cex.lab=2,axes=FALSE,xlab="",ylim=c(0,830))
lines(colMeans(sig[["clfswn"]],na.rm = TRUE),type="l",col="red",lwd=3)
axis(side = 1,labels = c("-2000","-1000","TSS","1000","2000"),at = c(1,13,26,39,50),lwd=2,cex=1.5,las=2,cex=2)
axis(side= 2, labels = seq(from=0,to=1000,by=200),at=seq(from=0,to=1000,by=200),lwd=2,cex=1.5)
mtext(text = "RPKM",side = 2,cex=2,line = 2.3,family="bold")
legend("topright",legend = c("Col-0",expression(italic("clfswn"))),col = c("blue","red"),lwd=3,cex=1)

## example plot for atbmi1abc vs Col-0
#plot(colMeans(sig[["Col-0"]],na.rm = TRUE),type="l",col="blue",lwd=3,ylab="",cex.lab=2,axes=FALSE,xlab="",ylim=c(0,200))
#lines(colMeans(sig[["atbmi1abc"]],na.rm = TRUE),type="l",col="red",lwd=3)
#axis(side = 1,labels = c("-2000","-1000","TSS","1000","2000"),at = c(1,13,26,39,50),lwd=2,cex=1.5,las=2,cex=2)
#axis(side= 2, labels = seq(from=0,to=200,by=40),at=seq(from=0,to=200,by=40),lwd=2,cex=1.5)
#mtext(text = "RPKM",side = 2,cex=2,line = 2.3,family="bold")
#legend("topright",legend = c("Col-0",expression(italic("atbmi1abc"))),col = c("blue","red"),lwd=3,cex=1)

## Storing the signal around tes
signal.around.tes <- sig

## Extracting signal for the entire gene
## Two different matrices will be used to store the signal for each gene
## in the two different conditions taken into account. (This will be slow)
## Expect an error at the end of the for loop and do not worry about it as it
## is controlled
final.result.col <- matrix(0,nrow=length(genes),ncol=100)
final.result.mut <- matrix(0,nrow=length(genes),ncol=100)

for(i in 1:length(genes))
{
  ## Printing current gene and extracting the features for the next two 
  ## genes
  print(i)
  current.gene.data <- genes.data[c(i,i+1),]
  
  ## Values for recentering (length of each feature divided by two)
  recenter.value <- floor(width(current.gene.data)/2)
  
  ## Centering features
  centered.current.gene.data <- current.gene.data
  start(centered.current.gene.data[1,]) <- start(centered.current.gene.data[1,]) + recenter.value[1]
  start(centered.current.gene.data[2,]) <- start(centered.current.gene.data[2,]) + recenter.value[2]
  end(centered.current.gene.data) <- start(centered.current.gene.data)
  
  start(current.gene.data[2,]) <- start(centered.current.gene.data[2,]) - recenter.value[1]
  end(current.gene.data[2,]) <- end(centered.current.gene.data[2,]) + recenter.value[1] -1
  
  ## Getting widths of the features
  widths <- width(current.gene.data)
  
  ## Getting the signal for the current gene data
  cvglists <- sapply(bigwig.files, import, 
                     format="BigWig", 
                     which=current.gene.data, 
                     as="RleList")
  ## You may want to change the names of the conditions here
  names(cvglists) <- c("Col-0", "clfswn")
  
  #names(cvglists) <- c("Col-0", "atbmi1abc") # example atbmi1abc vs Col-0
  if(widths[1] > 100)
  {
    sig <- featureAlignedSignal(cvglists, current.gene.data, #n.tile=100)
                                upstream=recenter.value[1], downstream=recenter.value[1],n.tile=100) 
    
    
    final.result.col[i,] <- sig[["Col-0"]][1,]
    final.result.mut[i,] <- sig[["clfswn"]][1,]
    
    #final.result.mut[i,] <- sig[["atbmi1abc"]][1,] # example atbmi1abc vs Col-0
  }
}

## Computing the average signal for each condition
mean.col <- colMeans(final.result.col)
mean.mut <- colMeans(final.result.mut)

## Merging the signal from the TSS, the gene body and TES
col.merged <- c(colMeans(signal.around.tss[["Col-0"]],na.rm = TRUE)[1:25],mean.col,colMeans(signal.around.tes[["Col-0"]],na.rm = TRUE)[26:50])
mut.merged <- c(colMeans(signal.around.tss[["clfswn"]],na.rm = TRUE)[1:25],mean.mut,colMeans(signal.around.tes[["clfswn"]],na.rm = TRUE)[26:50])

## mut.merged <- c(colMeans(signal.around.tss[["atbmi1abc"]],na.rm = TRUE)[1:25],mean.mut,colMeans(signal.around.tes[["atbmi1abc"]],na.rm = TRUE)[26:50])  ## example for atbmi1abc vs Col-0

## Saving the result in a file
data.col.mut <- data.frame(col.merged,mut.merged)
write.table(data.col.mut,row.names = F,col.names = F,file="metageneplot_data_clfswn_col0.txt")

## write.table(data.col.mut,row.names = F,col.names = F,file="metageneplot_data_atbmi1abc_col0.txt") ## example for atbmi1abc vs Col-0


## Plotting the final metageneplot
plot(col.merged,type="l",col="blue",lwd=3,ylab="",cex.lab=2.5,axes=FALSE,xlab="",ylim=c(0,1000))
lines(mut.merged,type="l",col="red",lwd=3)
axis(side = 1,labels = c("-2000","TSS","33%","66%","TES","2000"),at = c(1,25,58,91,125,150),lwd=2,cex=2,las=2,cex=2)
axis(side= 2, labels = seq(from=0,to=1000,by=200),at=seq(from=0,to=1000,by=200),lwd=2,cex=2)
mtext(text = "RPKM",side = 2,cex=2,line = 2.3,family="bold")
legend("topright",legend = c("Col-0",expression(italic("clfswn"))),col = c("blue","red"),lwd=3,cex=1)

## example atbmi1abc vs col-0
##plot(col.merged,type="l",col="blue",lwd=3,ylab="",cex.lab=2.5,axes=FALSE,xlab="",ylim=c(0,200))
##lines(mut.merged,type="l",col="red",lwd=3)
##axis(side = 1,labels = c("-2000","TSS","33%","66%","TES","2000"),at = c(1,25,58,91,125,150),lwd=2,cex=2,las=2,cex=2)
##axis(side= 2, labels = seq(from=0,to=200,by=40),at=seq(from=0,to=200,by=40),lwd=2,cex=2)
##mtext(text = "RPKM",side = 2,cex=2,line = 2.3,family="bold")
##legend("topright",legend = c("Col-0",expression(italic("atbmi1abc"))),col = c("blue","red"),lwd=3,cex=1)
