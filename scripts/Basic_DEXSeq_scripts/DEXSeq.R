# Rscript ./DEXSeq.R --count_matrix_dir PATH --gff_file_dir PATH --metadata PATH

args <- commandArgs( trailingOnly=TRUE )

if(!"DEXSeq" %in% rownames(installed.packages())){
  source("https://bioconductor.org/biocLite.R")
  biocLite("DEXSeq")
}

suppressPackageStartupMessages(library("DEXSeq"))
suppressPackageStartupMessages(library("dplyr"))

#load files
countFiles <- list.files(args[grep("--count_matrix_dir", args)+1], pattern = "fb.txt$", full.names = TRUE)
flattenedFile <- list.files(args[grep("--gff_file_dir", args)+1], pattern = "gff$", full.names = TRUE)
sampleTable <- read.table(args[grep("--metadata", args)+1], row.names = 1, header = TRUE, sep = ",")

#construct an DEXSeqDataSet object
dxd <- DEXSeqDataSetFromHTSeq(countfiles = countFiles, sampleData = sampleTable, design = ~ sample + exon + condition:exon, flattenedfile = flattenedfile)

BPPARAM = MultiCoreParam(workers=4)

# normalisation
dxd = estimateSizeFactors(dxd)

# dispertion estimation. estimate the variability of the data not explained by the biological variation.
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)

# test for each exon in each gene
dxd = testForDEU( dxd, BPPARAM=BPPARAM)

# calculate relative exon usage fold change
dxd = estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM)

# create a results object
dxr = DEXSeqResults( dxd )

dxr_dataframe = as.data.frame(dxr)

write.csv(dxr_dataframe, "DEE_results.csv")