#!/usr/bin/env Rscript
args = commandArgs(T)

print(dir())
#source("https://bioconductor.org/biocLite.R")
#biocLite('pasilla')
#biocLite("DESeq")

#library(pasilla)
#datafile = system.file( "extdata/pasilla_gene_counts.tsv", package="pasilla" )

library(DESeq)

CountTable = read.table(args[1], row.names = 1)

Design = read.table(args[2], row.names = 1)

cds = newCountDataSet(CountTable, Design)

cds = estimateSizeFactors( cds )

#sizeFactors( cds )
#counts( cds, normalized=TRUE )

cds = estimateDispersions( cds )
#View(cds)

plotDispEsts( cds )
#head( fData(cds) )
#View(cds)
res = nbinomTest( cds, "untreated", "treated" )
print(res)
plotMA(res)
