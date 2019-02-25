# Rscript EdgeR.R --condition string --counts PATH --metadata PATH
# code inspired by https://f1000research.com/articles/5-1408/v1

library(limma)
library(edgeR)
#library(dplyr)

args <- commandArgs( trailingOnly=TRUE )
args

if( "--condition" %in% args ){
  condition.idx <- grep("--condition", args)
  condition <- args[ condition.idx + 1 ]
} else {
  stop("please enter condition with prefix '--condition'")
}

if( "--counts" %in% args ){
  counts.idx <- grep("--counts", args)
  counts <- args[ counts.idx + 1 ]
  counts <- read.csv(counts, row.names = 1, header = TRUE)
} else {
  stop("please enter the path to the counts matrix with prefix '--counts'")
}

if( "--metadata" %in% args ){
  metadata.idx <- grep("--metadata", args)
  metadata <- args[ metadata.idx + 1 ]
  metadata <- read.csv(metadata, row.names = 1, header = TRUE)
} else {
  stop("please enter the path to the metadata with prefix '--metadata'")
}

DGE <- DGEList(counts = counts)
samplenames <- colnames(DGE)
samplenames

metadata <- metadata[samplenames,]
group <- as.factor(metadata[,condition])
DGE$samples$group <- group

#cpm <- cpm(DGE)

#keep.exprs <- rowSums(cpm>1)>=3
#DGE <- DGE[keep.exprs,, keep.lib.sizes=FALSE]

DGE <- calcNormFactors(DGE, method = "TMM")

design <- model.matrix(~0+group)
contr.matrix <- makeContrasts(contrasts = paste0("group0-group1"), levels = colnames(design))

v <- voom(DGE, design)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

dge.res <- topTable(efit,sort="none",n=Inf)
dge.res <- dge.res[,c("AveExpr", "logFC", "t", "P.Value", "adj.P.Val")]
colnames(dge.res) <- c("baseMean", "log2FoldChange", "stat", "pvalue","padj")
write.csv(dge.res, paste0("DGE_res.csv"))
