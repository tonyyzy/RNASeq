# Rscript featurecount.R --d PATH --g PATH --ft STRING --a STRING --mo LOGICAL --mm LOGICAL --f LOGICAL --lr LOGICAL --s INT --e LOGICAL --p INT

library("Rsubread")

args <- commandArgs( trailingOnly=TRUE )
args

if ("--d" %in% args){
  bam_dir.idx <- grep("--d", args)
  bam_dir.path <- args[ bam_dir.idx + 1 ]
  files <- list.files(bam_dir.path)
  files <- paste0(bam_dir.path,"/",files[grep(".bam", files)])
} else {
  stop("must provide directory with bam files as a parameter with prefix '--d'")
}

if ("--g" %in% args){
  gtf.idx <- grep("--g", args)
  gtf.path <- args[ gtf.idx + 1 ]
  gtf <- rtracklayer::import(gtf.path, format = "gtf")
  GTFAnnotationFile <- TRUE
} else {
  stop("please provide a gtf file with prefix '--g'")
}

if ("--ft" %in% args){
  featureType.idx <- grep("--ft", args)
  featureType <- args[ featureType.idx + 1 ]
} else {
  featureType <- "exon"
}

if ("--a" %in% args){
  attrType.idx <- grep("--a", args)
  attrType <- args[ attrType.idx + 1 ]
} else {
  attrType <- "gene_id"
}

if ("--mo" %in% args){
  MultipleOverlap <- TRUE
} else {
  MultipleOverlap <- FALSE
}

if ("--mm" %in% args){
  countMultiReads <- TRUE
} else {
  countMultiReads <- FALSE
}

if ("--f" %in% args){
  fraction <- TRUE
} else {
  fraction <- FALSE
}

if ("--lr" %in% args){
  LongReads <- TRUE
} else {
  LongReads <- FALSE
}

if ("--s" %in% args){
  stranded.idx <- grep("--s", args)
  stranded <- as.numeric(args[ stranded.idx + 1 ])
} else {
  stranded = 0
}

if ("--metadata" %in% args){
  metadata.idx <- grep("--e", args)
  metadata <- read.table(args[grep("--metadata", args)+1],header = TRUE, row.names = 1, sep=",")
  if("libType" %in% colnames(metadata)){
    libtype <- metadata$libType
  } else {
    stop("libtype must be a column in metadata")
  }
} else {
  stop("array of libtype must be provided with prefix '--e'")
}

if ("--p" %in% args){
  thread.idx <- grep("--p",args)
  thread <- as.numeric(args[ thread.idx + 1 ])
} else {
  thread <- 1
}

print("starting featurecounts")

for(i in 1:length(files)){
  x <- files[i]
  PairedEnd <- libtype[i]
  print(i)
  print(x)
  print(PairedEnd)
  if(PairedEnd == "paired-end"){
    PairedEnd <- TRUE
  } else if(PairedEnd == "single"){
    PairedEnd <- FALSE
  }
  if(exists("counts")){
    tmp_counts <- featureCounts(files=x, annot.ext = gtf.path, isGTFAnnotationFile = GTFAnnotationFile,
                          GTF.featureType = featureType, GTF.attrType = attrType,
                          allowMultiOverlap = MultipleOverlap, countMultiMappingReads = countMultiReads,
                          fraction = fraction, isLongRead = LongReads, strandSpecific = stranded,
                          isPairedEnd = PairedEnd, nthreads = thread)
    counts <- cbind(counts,tmp_counts$counts)
  } else {
    counts <- featureCounts(files=x, annot.ext = gtf.path, isGTFAnnotationFile = GTFAnnotationFile,
                         GTF.featureType = featureType, GTF.attrType = attrType,
                         allowMultiOverlap = MultipleOverlap, countMultiMappingReads = countMultiReads,
                         fraction = fraction, isLongRead = LongReads, strandSpecific = stranded,
                         isPairedEnd = PairedEnd, nthreads = thread)
    counts <- counts$counts
  }
}
colnames(counts) <- unlist(lapply(basename(files), function(x) gsub(".bam", "",x)))
write.csv(counts, "gene_count_matrix.csv")
