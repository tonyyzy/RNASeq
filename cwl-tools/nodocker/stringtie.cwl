#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: stringtie
inputs:
 input_bam:
  type: File
  inputBinding:
   position: 1
 threads:
   type: int
   inputBinding:
     position: 2
     prefix: -p
 annotation:
  type: File
  inputBinding:
   position: 3
   prefix: -G
 outfilename:
  type: string
  inputBinding:
    position: 4
    prefix: -o

outputs:
 example_out:
  type: File
  outputBinding:
   glob: $(inputs.outfilename)
