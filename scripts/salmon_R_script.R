# ./salmon_R_script.R --gtf FILE --metadata FILE --salmon_dir PATH

library(stringr)
suppressMessages(library(tximport))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(dplyr))
suppressMessages(library(rtracklayer))

args <- commandArgs( trailingOnly=TRUE )

if ("--gtf" %in% args) {
  gtf.file.idx  <- grep("--gtf", args)
  gtf.file.path <- args[ gtf.file.idx+1 ]
  if (file.exists(gtf.file.path)) {
    cat("\nLoading _serialised_ 'gtf' from \"",gtf.file.path,"\" ... ", sep="")
    TxDb <- makeTxDbFromGFF(gtf.file.path,format="gtf")
    cat("done\n")
  } else {
    cat("\nLocation of file containing 'gtf': \"",gtf.file.path,"\"\n", sep="")
    stop("File **DOES NOT EXIST**")
  }
}

if ("--metadata" %in% args) {
  metadata.file.idx  <- grep("--metadata", args)
  metadata.file.path <- args[ metadata.file.idx+1 ]
  if (file.exists(metadata.file.path)) {
    cat("\nLoading _serialised_ 'metadata' from \"",metadata.file.path,"\" ... ", sep="")
    samples <- read.table(metadata.file.path, header=TRUE, row.names = 1, sep = ",")
    cat("done\n")
  } else {
    cat("\nLocation of file containing 'metadata': \"",metadata.file.path,"\"\n", sep="")
    stop("File **DOES NOT EXIST**")
  }
}

if ("--salmon_dir" %in% args) {
  salmon_dir.idx  <- grep("--salmon_dir", args)
  salmon_dir.path <- args[ salmon_dir.idx+1 ]
  if (file.exists(salmon_dir.path)) {
    cat("\nLoading _serialised_ 'salmon_dir' from \"",salmon_dir.path,"\" ... ", sep="")
    files <- file.path(salmon_dir.path, rownames(samples), "quant.sf")
    cat("done\n")
  } else {
    cat("\nLocation of file containing 'salmon_dir': \"",salmon_dir.path,"\"\n", sep="")
    stop("File **DOES NOT EXIST**")
  }
}

if(nrow(samples) == length(files)){
  names(files) <- rownames(samples)
  k <- keys(TxDb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(TxDb, k,"GENEID", "TXNAME")
  tx2gene <- tx2gene[!is.na(tx2gene$GENEID),]

  txi <- tximport(files, type="salmon", tx2gene=tx2gene, dropInfReps=TRUE, ignoreTxVersion = TRUE)
  counts <- round(txi$counts,0)
  colnames(counts) <- names(files)
  write.csv(counts, "gene_count_matrix.csv")
  write.csv(txi$length, "gene_length_matrix.csv")
  write.csv(txi$abundance, "gene_abundance_matrix.csv")
} else {
  stop("Different number of samples, should be same number of salmon directories to samples in metadata")
}
