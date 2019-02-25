#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bioconductor-ballgown:2.14.0--r351_0

inputs:
  script:
    type: File
    inputBinding:
      position: 0
  tablemaker_output:
    type: Directory
    inputBinding:
      position: 1
      prefix: --data_dir
  metadata:
    type: File
    inputBinding:
      position: 2
      prefix: --metadata
  condition:
    type: string
    inputBinding:
      position: 3
      prefix: --condition

outputs:
  gene_matrix:
    type: File
    outputBinding:
      glob: "DGE_res.csv"
  transcript_matrix:
    type: File
    outputBinding:
      glob: "DTE_res.csv"
