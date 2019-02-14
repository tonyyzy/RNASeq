#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: tablemaker

hints:
  DockerRequirement:
    dockerPull: filipejesus/tablemaker:latest

arguments:
  - position: 2
    valueFrom: "-q"
  - position: 3
    valueFrom: "-W"

inputs:
  threads:
    type: string
    inputBinding:
      position: 1
      prefix: -q
  merged_gtf:
    type: File
    inputBinding:
      position: 3
      prefix: -G
  output:
    type: string
    inputBinding:
      position: 4
      prefix: -o
  bam:
    type: File
    inputBinding:
      position: 5

outputs:
  output:
    type: Directory
    inputBinding:
      glob: $(inputs.output)
