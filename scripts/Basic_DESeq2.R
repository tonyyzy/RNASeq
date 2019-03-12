# Rscript ./Basic_DESeq2.R --count_matrix PATH --metadata PATH --column name

#suppressMessages(library("dplyr"))

args <- commandArgs( trailingOnly=TRUE )

suppressMessages(library("DESeq2"))

data <- read.table(args[grep("--count_matrix", args)+1],header = TRUE, row.names = 1, sep=",", check.names = FALSE)
metadata <- read.table(args[grep("--metadata", args)+1],header = TRUE, row.names = 1, sep=",", check.names = FALSE)

if( "--condition" %in% args ){
  condition.idx <- grep("--condition", args)
  condition <- args[ condition.idx + 1 ]
  colnames(metadata) <- sub(condition, "condition",colnames(metadata))
}
print(metadata)
print(data)
comb <- combn(unique(as.character(metadata[,"condition"])), 2)
for(i in 1:ncol(comb)){
  metadata.f <- metadata[metadata$condition %in% comb[,i],]
  count.f <- data[,rownames(metadata.f)]
  dds <- DESeqDataSetFromMatrix(countData=count.f, colData=metadata.f, design= ~condition)

  dds <- estimateSizeFactors(dds)
  norm_count <- counts(dds, normalized=TRUE)

  dds <- DESeq(dds)
  
  res <- results(dds) 
  res <- data.frame(res@rownames,as.data.frame(res))
  colnames(res) <- c("name","basemean","log2foldchange","lfcSE","test_stat","p_value","p_adj")
  
  contrast <- gsub(".$","",paste0(paste0(unique(metadata.f$condition)),sep="-", collapse = ""))
  write.csv(res,paste0(contrast,"_","DGE_res.csv"), row.names = FALSE)
}

norm_count2 <- data.frame("name"=rownames(norm_count),as.data.frame(norm_count))
write.csv(norm_count2, "norm_count.csv", row.names = FALSE)
