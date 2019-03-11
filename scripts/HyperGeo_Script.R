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
    if(!"log2FoldChange" %in% colnames(Deseq.results)){
      stop("log2FoldChange must be a column in de results file")
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
Deseq.results <- Deseq.results[!is.na(Deseq.results$p_adj),]
list2 <- Deseq.results[Deseq.results$p_adj <= 0.05,]
results <- c()
for(list1 in Gene_Sets){
  overlap <- list1[list1 %in% rownames(list2)]
  tmp <- phyper(length(overlap),length(list1),(nrow(Deseq.results)-length(list1)),nrow(list2),lower.tail = FALSE, log.p = FALSE)
  results <- c(results,tmp)
}

print("step 4")
results2 <- p.adjust(results, method = "bonferroni")
results <- data.frame("pvalue"=results, "padj"=results2, row.names = names(Gene_Sets))
write.csv(results, "hypergeo_res.csv")
