#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: python2

inputs:
 program:
  type: File
  inputBinding:
    position: 1
 textSampleFile:
  type: File
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
