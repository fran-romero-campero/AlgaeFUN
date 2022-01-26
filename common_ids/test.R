description_table <- read.csv(file="ostreococcus_tauri_proteome_annotation.csv", header = F, sep="\t")
modified_table <- matrix(nrow=nrow(description_table), ncol=2)
modified_table[,1] <- description_table[,1]
i<-1
for( i in 1:nrow(description_table))
{
  modified_table[i,2] <- paste(description_table[i,2],
                               description_table[i,3],
                               description_table[i,4],
                               description_table[i,5],
                               description_table[i,6],
                               description_table[i,7],
                               description_table[i,8],
                               description_table[i,9],
                               description_table[i,10],
                               description_table[i,11],
                               description_table[i,12],
                               description_table[i,13],
                               description_table[i,14],
                               description_table[i,15],
                               description_table[i,16],
                               description_table[i,17],
                               description_table[i,18],
                               description_table[i,19],
                               description_table[i,20], sep=" ")
}
colnames(modified_table) <- c("Gene", "Common ID or description")
write.table(modified_table, file="otauri_IDs.tsv")
new_table<-read.table(file="otauri_IDs.tsv", header = T)
colnames(new_table) <- c("Genes", "Common ID or description")
target.genes <- c("ostta17g01460", "ostta02g00080", "ostta06g03840", "ostta04g01380")
conversion_d <- new_table$`Common ID or description`
names(conversion_d) <- new_table$Genes
final_conversion<-matrix(c(target.genes, conversion_d[target.genes]), ncol = 2, nrow = length(target.genes))
final.conversion <- as.data.frame(final_conversion)
colnames(final_conversion) <- c("Genes", "Common ID or description")

#primero fila despues columnas