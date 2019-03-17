# Rscript ./GSEA_Script.R --de_res PATH --gene_set PATH

args <- commandArgs( trailingOnly=TRUE )

if ("--gene_set" %in% args) {
  gene_set.file.idx  <- grep("--gene_set", args)
  gene_set.file.path <- args[ gene_set.file.idx+1 ]
  if (file.exists(gene_set.file.path)) {
    cat("\nLoading _serialised_ 'Gene Set' from \"",gene_set.file.path,"\" ... ", sep="")
    gs <- read.table(gene_set.file.path)
    gs <- gs[grep("HSA", gs$V2),]
    cat("done\n")
  } else {
    cat("\nLocation of file containing 'Gene Set': \"",gene_set.file.path,"\"\n", sep="")
    stop("File **DOES NOT EXIST**")
  }
}

## Read file containing the **serialised** 'SummarizedExperiment' object & check its class is correct ##
if (! "--de_res" %in% args) {
  cat("\n")
  stop("'Differential Expression' flag [--de_res] absent")
} else {
  de.file.path <- args[ grep("--de_res", args)+1 ]

  if (file.exists(de.file.path)) {
    cat("\nLoading _serialised_ 'differential expression' results from \"", de.file.path ,"\" ... ", sep="")
    Deseq.results <- read.csv(de.file.path, row.names = 1, header = TRUE)
    if(!"p_adj" %in% colnames(Deseq.results)){
      stop("p_adj must be a column in de results file")
    }
    if(!"log2foldchange" %in% colnames(Deseq.results)){
      stop("log2foldchange must be a column in de results file")
    }
    cat("done\n")
  }
  else {
    cat("\nLocation of file containing 'differential expression' results : \"", de.file.path,"\"\n", sep="")
    stop("File **DOES NOT EXIST**")
  }
}

print("step 1")

gs <- gs[gs[,1] %in% rownames(Deseq.results),]
Gene_Sets <- lapply(unique(gs[,2]), function(x) gs[gs[,2] == x,1])
names(Gene_Sets) <- unique(gs[,2])

print("step 2")

stats <- Deseq.results[!is.na(Deseq.results$p_adj),]
tmp <- stats[,"log2foldchange"]
names(tmp) <- rownames(stats)
print(head(tmp))

print("step 3")
gsea.Results <- fgsea::fgsea(Gene_Sets,tmp, 5000, minSize = 10)

print("step 4")

Results <- data.frame(gsea.Results[,2:7], row.names = gsea.Results$pathway)
write.csv(Results, "gsea_res.csv")
