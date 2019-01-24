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
  type: Directory
  inputBinding:
    position: 2
    prefix: -i

outputs:
  output:
    type: File[]
    outputBinding:
      glob: "*"

