#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript

requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bioconductor-ballgown:2.14.0--r351_0

arguments:
  - position: 4
    prefix: "&&"
    valueFrom: "ls"
    shellQuote: False

inputs:
  input_script:
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
  de_res:
    type: File[]
    outputBinding:
      glob: "*DGE_res.csv"
  ballgown_out:
    type: File[]
    outputBinding:
      glob: "*.csv"
