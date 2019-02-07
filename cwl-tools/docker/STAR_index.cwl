#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: STAR

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/star:2.6.0c--0

inputs:
  Threads:
    type: int
    inputBinding:
      prefix: --runThreadN
  Mode:
    type: string
    inputBinding:
      prefix: --runMode
  genomeDir:
    type: string
    default: "./"
    inputBinding:
      prefix: --genomeDir
  genomeFastaFiles:
    type: File
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
  genomeSAindexNbases:
    type: string?
    inputBinding:
      prefix: --genomeSAindexNbases
  RAMlimit:
    type: long
    default: 400000000000
    inputBinding:
      prefix: --limitGenomeGenerateRAM

outputs:
  output:
    type: File[]
    outputBinding:
      glob: "*"
