description_table <- read.csv(file="mpusilla_IDs.txt", header = F, sep="\t")
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
write.table(modified_table, file="mpusilla_IDs.tsv")
new_table<-read.csv(file="bprasinos_IDs.csv", header = T, sep="\t")

#primero fila despues columnas