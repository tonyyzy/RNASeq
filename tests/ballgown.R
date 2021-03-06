# Rscript ballgown.R --data_dir PATH --metadata PATH --condition STRING

library("matrixStats")
library("ballgown")

set.seed(1)

args <- commandArgs(trailingOnly = TRUE)

if( "--data_dir" %in% args){
  data_dir.idx <- grep("--data_dir", args)
  data_dir <- args[ data_dir.idx + 1 ]
} else {
  stop("please provide the path of the tablemaker output with prefix '--data_dir' \n")
}

if( "--metadata" %in% args ){
  metadata.idx <- grep("--metadata", args)
  metadata.path <- args[ metadata.idx + 1 ]
  if(file.exists(metadata.path)){
    metadata <- read.csv(metadata.path, header = TRUE)
  } else {
    stop("file <", metadata.path, "> not found. \n")
  }
} else {
  stop("please provide a metadata csv \n")
}

if( "--condition" %in% args ){
  condition.idx <- grep("--condition", args)
  condition <- args[ condition.idx + 1 ]
  if(!condition %in% colnames(metadata)){
    stop("condition given '", condition, "' not found in the metadata")
  }
} else {
  stop("please provide a valid condition with prefix '--conditions'")
}

comb <- combn(unique(as.character(metadata[,"condition"])), 2)
for(i in 1:ncol(comb)){
  metadata.f <- metadata[metadata$condition %in% comb[,i],]
  sample_full_path <- paste(data_dir, metadata.f[,1], sep = "/")
  bg <- ballgown(samples = as.vector(sample_full_path), pData = metadata.f)
  
  # Filter out transcripts with low variance
  bg_filt <- subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)
  
  # Perform DE
  results_transcripts <- stattest(bg_filt, feature="transcript", covariate=condition, getFC=TRUE, meas="FPKM")
  results_transcripts$fc <- log2(results_transcripts$fc)
  results_transcripts[,3:ncol(results_transcripts)] <- round(results_transcripts[,3:ncol(results_transcripts)],2)
  colnames(results_transcripts) <- c("feature", "name", "log2foldchange", "p_value","p_adj")
  results_transcripts <- results_transcripts[,c("name", "feature", "log2foldchange", "p_value","p_adj")]
  results_genes <- stattest(bg_filt, feature="gene", covariate=condition, getFC=TRUE, meas="FPKM")
  results_genes$fc <- log2(results_genes$fc)
  results_genes[,3:ncol(results_genes)] <- round(results_genes[,3:ncol(results_genes)],2)
  colnames(results_genes) <- c("feature", "name", "log2foldchange", "p_value","p_adj")
  results_genes <- results_genes[,c("name", "feature", "log2foldchange", "p_value","p_adj")]
  contrast <- gsub(".$","",paste0(paste0(unique(metadata.f$condition)),sep="-", collapse = ""))
  write.csv(results_transcripts,paste0(contrast,"_","DTE_res.csv"), row.names = FALSE)
  write.csv(results_genes,paste0(contrast,"_","DGE_res.csv"), row.names = FALSE)
  norm_count <- gexpr(bg)
  norm_count <- data.frame("name"=rownames(norm_count), norm_count)
  write.csv(norm_count,paste0(contrast,"_norm_count.csv"), row.names = FALSE)
}
