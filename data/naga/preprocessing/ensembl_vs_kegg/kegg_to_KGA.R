naga.xp <- read.table(file="result_identity.txt",as.is=T)
head(naga.xp)

naga.info <- read.table(file="GCF_000240725.1_ASM24072v1_feature_table.txt",header=T,sep="\t",as.is=T)
head(naga.info)

naga.info.filtered <- subset(naga.info, feature == "mRNA")

head(naga.info.filtered)
xp.ids <- naga.info.filtered$related_accession
nga.ids <- naga.info.filtered$locus_tag

names(nga.ids) <- xp.ids
nga.ids[1:3]

naga.xp.res <- read.table(file="result_identity.txt",as.is=T,header=F)

naga.ids <- naga.xp.res$V1
nga.xp <- nga.ids[naga.xp.res$V2]
names(nga.xp) <- NULL
sum(is.na(nga.xp))

naga.nga.correspondence <- data.frame(naga=naga.ids, nga=nga.xp,stringsAsFactors = F)

write.table(x = naga.nga.correspondence,file = "naga_nga_correspondence.txt",quote = F,sep = "\t",row.names = F,col.names = F)
