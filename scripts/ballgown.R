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
sample_full_path <- paste(data_dir, metadata[,1], sep = "/")
bg <- ballgown(samples = as.vector(sample_full_path), pData = metadata)

# Filter out transcripts with low variance
bg_filt <- subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Perform DE
results_transcripts <- stattest(bg_filt, feature="transcript", covariate=condition, getFC=TRUE, meas="FPKM")
results_transcripts$fc <- log2(results_transcripts$fc)
colnames(results_transcripts) <- c("feature", "name", "log2foldchage", "p_value","p_adj")
results_genes <- stattest(bg_filt, feature="gene", covariate=condition, getFC=TRUE, meas="FPKM")
results_genes$fc <- log2(results_genes$fc)
colnames(results_genes) <- c("feature", "name", "log2foldchage", "p_value","p_adj")
write.csv(results_transcripts,"DTE_res.csv", row.names = FALSE)
write.csv(results_genes,"DGE_res.csv", row.names = FALSE)
