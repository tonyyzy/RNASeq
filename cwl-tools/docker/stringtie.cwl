#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: stringtie
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/stringtie:1.3.0--hd28b015_2
inputs:
 input_bam:
  type: File
  inputBinding:
   position: 1
   prefix: -eB
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
 stringtie_out:
  type: File
  outputBinding:
   glob: $(inputs.outfilename)
