#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: $(inputs.bam_files)
  DockerRequirement:
    dockerPull: filipejesus/featurecounts:3.8

arguments:
  - prefix: --d
    position: 1
    valueFrom: $(runtime.outdir)

inputs:
  input_script:
    type: File
    inputBinding:
      position: 0
  bam_files:
    type: File[]
  gtf:
    type: File
    inputBinding:
      position: 2
      prefix: --g
  threads:
    type: int?
    inputBinding:
      position: 3
      prefix: --p
  featuretype:
    type: string?
    inputBinding:
      position: 4
      prefix: --ft
  attribute:
    type: string?
    inputBinding:
      position: 5
      prefix: --a
  multipleoverlap:
    type: string?
    inputBinding:
      position: 6
      prefix: --mo
  multiread:
    type: string?
    inputBinding:
      position: 7
      prefix: --mm
  fraction:
    type: string?
    inputBinding:
      position: 8
      prefix: --f
  longreads:
    type: string?
    inputBinding:
      position: 9
      prefix: --lr
  stranded:
    type: string?
    inputBinding:
      position: 10
      prefix: --s
  metadata:
    type: File
    inputBinding:
      position: 11
      prefix: --metadata

outputs:
  gene_count_output:
    type: File
    outputBinding:
      glob: "gene_count_matrix.csv"
