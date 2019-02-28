#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
hints:
 DockerRequirement:
  dockerPull: biowardrobe2/scidap-deseq:v0.0.5
requirements:
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
    glob: "*DGE_results.csv"
