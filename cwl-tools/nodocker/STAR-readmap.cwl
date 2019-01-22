#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: STAR

inputs:
  Threads:
    type: string
    inputBinding:
      prefix: --runThreadN
  genomeDir:
    type: Directory
    inputBinding:
      prefix: --genomeDir
  readFilesIn:
    type: File[]?
    inputBinding:
      prefix: --readFilesIn
  outFileNamePrefix:
    type: string?
    inputBinding:
      prefix: --outFileNamePrefix

outputs:
  output:
    type: File[]
    outputBinding:
      glob: "*"
    
    