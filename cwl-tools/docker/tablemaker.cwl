#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: tablemaker

hints:
  DockerRequirement:
    dockerPull: filipejesus/cufflinks:latest

requirements:
   InlineJavascriptRequirement: {}
   ShellCommandRequirement: {}

arguments:
  - position: 2
    valueFrom: "-q"
  - position: 3
    valueFrom: "-W"

inputs:
  threads:
    type: int
    inputBinding:
      position: 1
      prefix: -p
  merged_gtf:
    type: File
    inputBinding:
      position: 4
      prefix: -G
  output:
    type: string
    inputBinding:
      position: 5
      prefix: -o
  bam:
    type: File
    inputBinding:
      position: 6

outputs:
  tablemaker_out:
    type: Directory
    outputBinding:
      glob: $(inputs.output)
