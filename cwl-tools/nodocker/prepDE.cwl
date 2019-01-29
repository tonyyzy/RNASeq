#!/usr/bin/env cwl-runner


cwlVersion: v1.0
class: CommandLineTool
baseCommand: python2
inputs:
 program:
  type: string
  inputBinding:
    position: 1
 gtfDir:
  type: Directory[]
  inputBinding:
    position: 2
    prefix: -i

outputs:
  gene_output:
    type: File
    outputBinding:
      glob: "*gene_count_matrix*"
  transcript_output:
    type: File
    outputBinding:
      glob: "*transcript_count_matrix*"
