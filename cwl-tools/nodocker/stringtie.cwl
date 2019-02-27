#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: stringtie
inputs:
 bam:
  type: File
  inputBinding:
   position: 1
   prefix: -eB
 threads:
   type: int
   inputBinding:
     position: 2
     prefix: -p
 gtf:
  type: File
  inputBinding:
   position: 3
   prefix: -G
 output:
  type: string
  inputBinding:
    position: 4
    prefix: -o

outputs:
 stringtie_out:
  type: File
  outputBinding:
   glob: $(inputs.output)
