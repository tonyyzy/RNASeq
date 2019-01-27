#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: STAR

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/star:2.6.0c--0

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
    type: File[]
    inputBinding:
      prefix: --readFilesIn
  outFileNamePrefix:
    type: string?
    inputBinding:
      prefix: --outFileNamePrefix

outputs:
  sam_output:
    type: File
    outputBinding:
      glob: "*.sam"
  log_outputs:
    type: File[]
    outputBinding:
      glob: "*.out"
  tab_output:
    type: File
    outputBinding:
      glob: "*.tab"
    
    
