#/usr/bin/env Rscript

test_path <- system.file("tests", package = "rtracklayer")
test_gff3 <- file.path(test_path, "genes.gff3")
test <- import(test_gff3)
export(test,"test.gtf","gtf")