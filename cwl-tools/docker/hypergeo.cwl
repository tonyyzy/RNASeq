#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript

requirements:
 DockerRequirement:
  dockerPull: biowardrobe2/scidap-deseq:v0.0.5
  
inputs:
 input_script:
  type: File
  inputBinding:
    position: 1
 de_res:
  type: File
  inputBinding:
   position: 2
   prefix: --de_res
 gene_set:
   type: File
   inputBinding:
     position: 3
     prefix: --gene_set

outputs:
  hypergeo_out:
   type: File
   outputBinding:
    glob: "hypergeo_res.csv"
