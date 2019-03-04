#/usr/bin/env RScript

# Rscript ./DEXSeq.R --count_matrix_dir PATH --gff_file_dir PATH --metadata PATH

args <- commandArgs( trailingOnly=TRUE )

suppressPackageStartupMessages(library("DEXSeq"))

#load files
countFiles <- list.files(args[grep("--count_matrix_dir", args)+1], pattern = "htseq_count.csv$", full.names = TRUE)
flattenedFile <- list.files(args[grep("--gff_file_dir", args)+1], pattern = "gff$", full.names = TRUE)
sampleTable <- read.table(args[grep("--metadata", args)+1], row.names = 1, header = TRUE, sep = ",")

if( "--condition" %in% args ){
  condition.idx <- grep("--condition", args)
  condition <- args[ condition.idx + 1 ]
  colnames(sampleTable) <- sub(condition, "condition",colnames(sampleTable))
}

#construct an DEXSeqDataSet object
dxd <- DEXSeqDataSetFromHTSeq(countfiles = countFiles, sampleData = sampleTable, design = ~ sample + exon + condition:exon, flattenedfile = flattenedFile)

if("--threads" %in% args){
  threads.idx <- grep("--threads", args)
  threads <- args[ threads.idx + 1 ]
  # Set up workers
  BPPARAM = MultiCoreParam(workers = as.numeric(threads))
  # normalisation
  dxd = estimateSizeFactors( dxd )
  # dispertion estimation. estimate the variability of the data not explained by the biological variation.
  dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)
  # test for each exon in each gene
  dxd = testForDEU( dxd, BPPARAM=BPPARAM)
  # calculate relative exon usage fold change
  dxd = estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)
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

dxd_norm <- counts(dxd,normalized=T)

# create a results object
dxr = DEXSeqResults( dxd )

dxr_dataframe = as.data.frame(dxr)

write.csv(dxr_dataframe, "DEE_results.csv")
dxd_norm2 <- data.frame("name"=rownames(dxd_norm),as.data.frame(dxd_norm))
write.csv(dxd_norm2, "norm_count.csv", row.names = FALSE)
