#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand:
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/star:2.6.0c--0
arguments:
  - position: 1
    valueFrom: "mkdir"
  - position: 2
    valueFrom: $(inputs.output)
  - position: 3
    valueFrom: "&&"
    shellQuote: False
  - position: 4
    valueFrom: "STAR"

inputs:
  threads:
    type: int
    inputBinding:
      position: 5
      prefix: --runThreadN
  Mode:
    type: string
    inputBinding:
      position: 6
      prefix: --runMode
  output:
    type: string
    default: "STARIndex"
    inputBinding:
      position: 7
      prefix: --genomeDir
  fasta:
    type: File[]
    inputBinding:
      position: 8
      prefix: --genomeFastaFiles
  gtf:
    type: File[]
    inputBinding:
      position: 9
      prefix: --sjdbGTFfile
  sjdbGTFtagExonParentTranscript:
    type: string?
    inputBinding:
      position: 10
      prefix: --sjdbGTFtagExonParentTranscript
  genomeSAindexNbases:
    type: string?
    inputBinding:
      position: 11
      prefix: --genomeSAindexNbases
  RAMlimit:
    type: long
    default: 400000000000
    inputBinding:
      position: 12
      prefix: --limitGenomeGenerateRAM

outputs:
  index_out:
    type: Directory
    outputBinding:
      glob: $(inputs.output)
