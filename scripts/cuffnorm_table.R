args <- commandArgs( trailingOnly=TRUE )

metadata <- read.csv(args[grep("--metadata", args)+1], row.names = 1, check.names = FALSE, header = TRUE)
count_table <- read.table(args[grep("--norm_table", args)+1], sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
gene_info <- read.table(args[grep("--attr_table", args)+1], sep = "\t", header = TRUE, check.names = FALSE)

count_table <- count_table[gene_info$tracking_id,]
rownames(count_table) <- gene_info$gene_short_name

columns <- c()
for(x in unique(metadata$condition)){
  tmp.metadata <- metadata[metadata$condition == x,]
  columns <- c(columns, rownames(tmp.metadata))
}

colnames(count_table) <- columns
counts <- data.frame("name" = rownames(count_table), count_table)
write.csv(counts, paste0(args[grep("--output", args)+1],"/norm_count.csv"), row.names = FALSE)