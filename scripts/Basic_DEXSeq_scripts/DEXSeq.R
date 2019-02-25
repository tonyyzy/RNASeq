#/usr/bin/env RScript

# Rscript ./DEXSeq.R --count_matrix_dir PATH --gff_file_dir PATH --metadata PATH

args <- commandArgs( trailingOnly=TRUE )

suppressPackageStartupMessages(library("DEXSeq"))

#load files
countFiles <- list.files(args[grep("--count_matrix_dir", args)+1], pattern = "htseq_count.csv$", full.names = TRUE)
flattenedFile <- list.files(args[grep("--gff_file_dir", args)+1], pattern = "gff$", full.names = TRUE)
sampleTable <- read.table(args[grep("--metadata", args)+1], row.names = 1, header = TRUE, sep = ",")

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

# create a results object
dxr = DEXSeqResults( dxd )

dxr_dataframe = as.data.frame(dxr)

write.csv(dxr_dataframe, "DEE_results.csv")
