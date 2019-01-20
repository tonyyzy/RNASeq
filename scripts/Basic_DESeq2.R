# Rscript ./Basic_DESeq2.R --count_matrix PATH --metadata PATH --column name

suppressMessages(library("dplyr"))

args <- commandArgs( trailingOnly=TRUE )

if(installed.packages() %>% .[,1] %>% any(grepl("DESeq2",.)) == FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")  
}

suppressMessages(library("DESeq2"))

data <- read.table(args[grep("--count_matrix", args)+1],header = TRUE, row.names = 1, sep=",")
metadata <- read.table(args[grep("--metadata", args)+1],header = TRUE, row.names = 1, sep=",")

metadata <- metadata[colnames(data),]

dds <- DESeqDataSetFromMatrix(countData=data, colData=metadata, design= ~condition) %>% DESeq(.)

res <- results(dds) %>% as.data.frame(.)

write.csv(res,"DGE_results.csv")
