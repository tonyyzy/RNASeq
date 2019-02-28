# Rscript ./Basic_DESeq2.R --count_matrix PATH --metadata PATH --column name

#suppressMessages(library("dplyr"))

args <- commandArgs( trailingOnly=TRUE )

suppressMessages(library("DESeq2"))

data <- read.table(args[grep("--count_matrix", args)+1],header = TRUE, row.names = 1, sep=",")
metadata <- read.table(args[grep("--metadata", args)+1],header = TRUE, row.names = 1, sep=",")

metadata <- metadata[colnames(data),]

if( "--condition" %in% args ){
  condition.idx <- grep("--condition", args)
  condition <- args[ condition.idx + 1 ]
  colnames(metadata) <- sub(condition, "condition",colnames(metadata))
}

comb <- combn(unique(metadata[,"condition"]), 2)
for(i in 1:ncol(comb)){
  metadata.f <- metadata[metadata$condition %in% comb[,i],]
  count.f <- data[,rownames(metadata.f)]
  dds <- DESeqDataSetFromMatrix(countData=count.f, colData=metadata.f, design= ~condition)
  dds <- DESeq(dds)
  
  res <- results(dds) 
  res <- as.data.frame(res)
  
  contrast <- gsub(".$","",paste0(paste0("group", unique(metadata$condition)),sep="-", collapse = ""))
  write.csv(res,paste0(contrast,"_","DGE_results.csv"))
}
