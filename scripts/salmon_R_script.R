# ./salmon_R_script.R --gtf FILE --metadata FILE --salmon_dir PATH

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
  exon_version_df <- data.frame("WithVersion"=NA,"TXNAME"=NA)
  for(x in 1:length(files)){
    sf_table <- read.table(files[x], header=TRUE)
    exon_df <- data.frame("WithVersion"=sf_table$Name,"TXNAME"=gsub("..$", "", sf_table$Name))
    exon_version_df <- rbind(exon_version_df,exon_df)
    exon_version_df <- exon_version_df[!duplicated(exon_version_df$WithVersion),]
  }
  exon_version_df <- exon_version_df[!is.na(exon_version_df$WithVersion),]
  
  k <- keys(TxDb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(TxDb, k,"GENEID", "TXNAME")
  tx2gene <- left_join(exon_version_df,tx2gene, by = "TXNAME")
  tx2gene <- tx2gene[!is.na(tx2gene$GENEID),]
  tx2gene$TXNAME <- NULL
  colnames(tx2gene) <- c("TXNAME", "GENEID")
  
  txi <- tximport(files, type="salmon", tx2gene=tx2gene, dropInfReps=TRUE)
  write.csv(txi$counts, "gene_count_matrix.csv")
} else {
  stop("Different number of samples, should be same number of salmon directories to samples in metadata")
}
