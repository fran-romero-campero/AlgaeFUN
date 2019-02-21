BiocManager::install("clusterProfiler", version = "3.8")
library(clusterProfiler)

library(org.Otauri.eg.db)

ostta.example <- read.table(file = "ostta_example_2_activated_genes.txt",header = FALSE,as.is = TRUE)[[1]]
length(ostta.example)

ostta.universe <- unique(select(org.Otauri.eg.db,columns = c("GO"),keys=keys(org.Otauri.eg.db,keytype = "GID"))[["GID"]])
length(ostta.universe)

ego <- enrichGO(gene          = ostta.example,
                universe      = ostta.universe,
                OrgDb         = org.Otauri.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE,
                keyType = "GID")


barplot(ego,drop=TRUE,showCategory = 10)

res.ego <- as.data.frame(ego)
head(res.ego)

png(filename = "goplot.png",width = 1000,height = 1000,res = 100)
goplot(ego,showCategory = 10)
dev.off()

dotplot(ego)

emapplot(ego)

cnetplot(ego)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview", version = "3.8")
library("pathview")

kk <- enrichKEGG(gene = paste0("OT_",ostta.example), organism = "ota",keyType = "kegg",
                 universe = paste0("OT_",ostta.universe),qvalueCutoff = 0.05)

head(kk)
res.kk <- as.data.frame(kk)

mkk <- enrichMKEGG(gene = paste0("OT_",ostta.example), organism = "ota",keyType = "kegg")
head(mkk)

genes.pathway <- rep(0,length(ostta.universe))
names(genes.pathway) <- paste0("OT_",ostta.universe)

genes.pathway[paste0("OT_",ostta.example)] <- 1

for(i in 1:nrow(res.kk))
{
  pathview(gene.data = sort(genes.pathway,decreasing = TRUE), pathway.id = res.kk$ID[i], species = "ota",limit = list(gene=max(abs(genes.pathway)), cpd=1),gene.idtype ="kegg")
}

pathview(gene.data = sort(genes.pathway,decreasing = TRUE), pathway.id = "ota03030", species = "ota",limit = list(gene=max(abs(genes.pathway)), cpd=1),gene.idtype ="kegg")
