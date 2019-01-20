#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand:
hints:
  DockerRequirement:
    dockerPull: dceoy/star
requirements:
  ResourceRequirement:
    tmpdirMin: 100000
    outdirMin: 100000

inputs:
  Threads:
    type: string
    inputBinding:
      position: 1
      prefix: --runThreadN
  Mode:
    type: string
    inputBinding:
      position: 2
      prefix: --runMode
  genomeDir:
    type: string
    inputBinding:
      position: 3
      prefix: --genomeDir
  genomeFastaFiles:
    type: File
    inputBinding:
      position: 4
      prefix: --genomeFastaFiles
  sjdbGTFfile:
    type: File
    inputBinding:
      position: 5
      prefix: --sjdbGTFfile
  sjdbGTFtagExonParentTranscript:
    type: string
    inputBinding:
      position: 6
      prefix: --sjdbGTFtagExonParentTranscript
  sjdbOverhang:
    type: string
    inputBinding:
      position: 7
      prefix: --sjdbOverhang
  
outputs:
  output:
    type: File[]
    outputBinding:
      glob: "*"
