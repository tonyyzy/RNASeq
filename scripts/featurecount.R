# Rscript featurecount.R --d PATH --g PATH --ft STRING --a STRING --mo LOGICAL --mm LOGICAL --f LOGICAL --lr LOGICAL --s INT --e LOGICAL --p INT

library("Rsubread")

args <- commandArgs( trailingOnly=TRUE )
args

if ("--d" %in% args){
  bam_dir.idx <- grep("--d", args)
  bam_dir.path <- args[ bam_dir.idx + 1 ]
  files <- list.files(bam_dir.path)
} else {
  stop("must provide directory with bam files as a parameter with prefix '--d'")
}

if ("--g" %in% args){
  gtf.idx <- grep("--g", args)
  gtf <- args[ gtf.idx + 1 ]
  gtf <- rtracklayer::import(gtf, format = "gtf")
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

if ("--e" %in% args){
  PairedEnd <- TRUE
} else {
  PairedEnd <- FALSE
}

if ("--p" %in% args){
  thread.idx <- grep("--p",args)
  thread <- as.numeric(args[ thread.idx + 1 ])
} else {
  thread <- 1
}

for(x in files){
  tmp <- featureCounts(files=x, annot.ext = gtf, isGTFAnnotationFile = GTFAnnotationFile,
                       GTF.featureType = featureType, GTF.attrType = attrType, 
                       allowMultiOverlap = MultipleOverlap, countMultiMappingReads = countMultiRead, 
                       fraction = fraction, isLongRead = LongRead, strandSpecific = stranded, 
                       isPairedEnd = PairedEnd, nthreads = thread)
}
