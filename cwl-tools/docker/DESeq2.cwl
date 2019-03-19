#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
requirements:
 DockerRequirement:
  dockerPull: biowardrobe2/scidap-deseq:v0.0.5
 InlineJavascriptRequirement: {}

baseCommand: Rscript
inputs:
 input_script:
  type: File
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
 threads:
   type: int?
   inputBinding:
     position: 4
     prefix: --threads
outputs:
  DESeq2_out:
   type: File[]
   outputBinding:
    glob: "*.csv"
  de_res:
    type: File[]
    outputBinding:
      glob: "*DGE_res.csv"
