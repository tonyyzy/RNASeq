#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript
inputs:
 script:
  type: string
  default: "/project/home18/apg3718/RNASeq/scripts/Basic_DESeq2.R"
  inputBinding:
    position: 1
 count_matrix:
  type: File
  inputBinding:
   position: 2
   prefix: --count_matrix
 metadata:
   type: File
   inputBinding:
     position: 3
     prefix: --metadata
outputs:
  DESeq2_out:
   type: File
   outputBinding:
    glob: DGE_results.csv
