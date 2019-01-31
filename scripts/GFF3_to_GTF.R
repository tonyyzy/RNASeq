#/usr/bin/env Rscript --gff3 PATH --outname PATH

args <- commandArgs( trailingOnly=TRUE )

library(rtracklayer)

gff3.file.path <- args[ grep("--gff3", args)+1 ]
output.name <- args[ grep("--outname", args)+1 ]

test <- import(con = gff3.file.path)
export.gff2(object = test,con = output.name)