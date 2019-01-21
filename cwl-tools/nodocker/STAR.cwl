#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: STAR

inputs:
  Threads:
    type: string
    inputBinding:
      prefix: --runThreadN
  Mode:
    type: string
    inputBinding:
      prefix: --runMode
  genomeDir:
    type: string
    inputBinding:
      prefix: --genomeDir
  genomeFastaFiles:
    type: File?
    inputBinding:
      prefix: --genomeFastaFiles
  sjdbGTFfile:
    type: File
    inputBinding:
      prefix: --sjdbGTFfile
  sjdbGTFtagExonParentTranscript:
    type: string?
    inputBinding:
      prefix: --sjdbGTFtagExonParentTranscript
  sjdbOverhang:
    type: string?
    inputBinding:
      prefix: --sjdbOverhang
  genomeSAindexNbases:
    type: string?
    inputBinding:
      prefix: --genomeSAindexNbases
  readFilesIn:
    type: string?
    inputBinding:
      prefix: --readFilesIn
  outSAMtype:
    type: string?
    inputBinding:
      prefix: --outSAMtype
  outFileNamePrefix:
    type: string?
    inputBinding:
      prefix: --outFileNamePrefix

outputs:
  output:
    type: File[]
    outputBinding:
      glob: "*"
