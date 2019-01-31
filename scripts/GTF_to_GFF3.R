#/usr/bin/env Rscript --gtf PATH --outname PATH
args <- commandArgs( trailingOnly=TRUE )

library(rtracklayer)

gtf.file.path <- args[ grep("--gtf", args)+1 ]
output.name <- args[ grep("--outname", args)+1 ]

test <- import(gtf.file.path, "gtf")
export(test,output.name,"gff3")