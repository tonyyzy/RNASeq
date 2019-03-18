#/usr/bin/env RScript

# Rscript ./DEXSeq.R --count_matrix_dir PATH --gff_file_dir PATH --metadata PATH

args <- commandArgs( trailingOnly=TRUE )

suppressPackageStartupMessages(library("DEXSeq"))
suppressPackageStartupMessages(library("BiocParallel"))
suppressPackageStartupMessages(library("dplyr"))

set.seed(1)

#load files
countFiles <- list.files(args[grep("--count_matrix_dir", args)+1], pattern = "htseq_count.csv$", full.names = TRUE)
metadata <- read.table(args[grep("--metadata", args)+1], row.names = 1, header = TRUE, sep = ",")

if( "--condition" %in% args ){
  condition.idx <- grep("--condition", args)
  condition <- args[ condition.idx + 1 ]
  colnames(metadata) <- sub(condition, "condition",colnames(metadata))
}

df <- read.table(countFiles[1], sep = "\t")
df <- df[!df$V1 %in% c("_ambiguous", "_ambiguous_readpair_position", "_empty", "_lowaqual", "_notaligned"),]
colnames(df) <- c("V1", sub("_htseq_count.csv","",basename(countFiles[1])))
for (x in countFiles[2:length(countFiles)]){
  tmp.df <- read.table(x, sep = "\t")
  tmp.df <- tmp.df[!tmp.df$V1 %in% c("_ambiguous", "_ambiguous_readpair_position", "_empty", "_lowaqual", "_notaligned"),]
  colnames(tmp.df) <- c("V1", sub("_htseq_count.csv","",basename(x)))
  df <- full_join(df, tmp.df, by = "V1")
}
df[is.na(df)] <- 0
rownames(df) <- df$V1
df$V1 <- NULL

comb <- combn(unique(as.character(metadata[,"condition"])), 2)
for(i in 1:ncol(comb)){
  metadata.f <- metadata[metadata$condition %in% comb[,i],]
  count.f <- df[,rownames(metadata.f)]
  #construct an DEXSeqDataSet object
  ids <- t(as.data.frame(strsplit(rownames(count.f),":")))
  dxd <- DEXSeqDataSet(count.f, metadata, design= ~ sample + exon + condition:exon, featureID = ids[,2], groupID = ids[,1])

  if("--threads" %in% args){
    threads.idx <- grep("--threads", args)
    threads <- args[ threads.idx + 1 ]
    # Set up workers
    BPPARAM = BiocParallel::MulticoreParam(workers = as.numeric(threads))
    # normalisation
    dxd = estimateSizeFactors( dxd )
    # dispertion estimation. estimate the variability of the data not explained by the biological variation.
    dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)
    # test for each exon in each gene
    dxd = testForDEU( dxd, BPPARAM=BPPARAM)
    # calculate relative exon usage fold change
    dxd = estimateExonFoldChanges(dxd, BPPARAM=BPPARAM, fitExpToVar="condition")
  } else {
    # normalisation
    dxd = estimateSizeFactors(dxd)
    # dispertion estimation. estimate the variability of the data not explained by the biological variation.
    dxd = estimateDispersions( dxd)
    # test for each exon in each gene
    dxd = testForDEU( dxd)
    # calculate relative exon usage fold change
    dxd = estimateExonFoldChanges(dxd, fitExpToVar="condition")
  }
  # create a results object
  contrast <- gsub(".$","",paste0(paste0(unique(metadata.f$condition)),sep="-", collapse = ""))

  dxr = DEXSeqResults( dxd )
  dxr_dataframe = as.data.frame(dxr)
  dxr_dataframe[,c(3,4,5,6,7,8,9,10)]<-round(dxr_dataframe[,c(3,4,5,6,7,8,9,10)],0)
  dxr_dataframe <- data.frame("name"=rownames(dxr_dataframe),dxr_dataframe)
  dxr_dataframe <- dxr_dataframe[,1:12]
  colnames(dxr_dataframe) <- c("name", "groupID", "featureID","exonBaseMean","dispersion","test_stat","p_value","p_adj","Normal","Tumour","log2foldchange","genomicData")
  write.csv(dxr_dataframe, paste0(contrast,"_DEE_results.csv"), row.names = FALSE)

  dxd_norm <- counts(dxd,normalized=T)
  dxd_norm2 <- data.frame("name"=rownames(dxd_norm),as.data.frame(dxd_norm))
  write.csv(dxd_norm2, paste0(contrast,"_norm_count.csv"), row.names = FALSE)
}
