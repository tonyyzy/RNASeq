library(tximport)
library(GenomicAlignments)
library(GenomicFeatures)
library(dplyr)

TxDb <- makeTxDbFromGFF(file,format=c("auto", "gff3", "gtf"))

samples <- read.table(metadata, header=TRUE, row.names = 1, sep = ",")
files <- file.path(salmon_dir, rownames(samples), "quant.sf")
names(files) <- rownames(samples)
tmp4 <- data.frame("WithVersion"=NA,"TXNAME"=NA)
for(x in 1:length(files)){
  tmp2 <- read.table(files[x], header=TRUE)
  tmp3 <- data.frame("WithVersion"=tmp2$Name,"TXNAME"=gsub("..$", "", tmp2$Name))
  tmp4 <- rbind(tmp4,tmp3)
  tmp4 <- tmp4[!duplicated(tmp4$WithVersion),]
}
tmp4 <- tmp4[!is.na(tmp4$WithVersion),]

TxDb <- makeTxDbFromGFF(gtf_file,format=c("gtf"))
k <- keys(TxDb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(TxDb, k,"GENEID", "TXNAME")
tx2gene <- left_join(tmp4,tx2gene, by = "TXNAME")
tx2gene <- tx2gene[!is.na(tx2gene$GENEID),]
tx2gene$TXNAME <- NULL
colnames(tx2gene) <- c("TXNAME", "GENEID")

txi <- tximport(files, type="salmon", tx2gene=tx2gene, dropInfReps=TRUE)
write.csv(txi$counts, paste0(salmon_dir,"/count_matrix.csv"))