#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript

hints:
  DockerRequirement:
    dockerPull: filipejesus/featurecounts:latest

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: $(inputs.bam_files)

arguments:
  - prefix: --d
    position: 1
    valueFrom: $(runtime.outdir)

inputs:
  script:
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
    type: string?
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
  pairedend:
    type: string?
    inputBinding:
      position: 11
      prefix: --e

outputs:
  output:
    type: File
    outputBinding:
      glob: $("gene_count_matrix.csv")
