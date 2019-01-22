#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
requirements:
  - class: ShellCommandRequirement
hints:
  DockerRequirement:
    dockerPull: genomicpariscentre/samtools
baseCommand: samtools
inputs:
 action:
  type: string
  default: "view"
  inputBinding:
   position: 1
   shellQuote: False
 threads:
  type: int
  inputBinding:
    prefix: -@
    position: 2
    shellQuote: False
 samfile:
   type: File
   inputBinding:
    position: 3
    prefix: -Su
    shellQuote: False
 pipe:
    type: string
    default: "|"
    inputBinding:
      position: 4
      shellQuote: False
 samtools:
  type: string
  default: "samtools"
  inputBinding:
   position: 5
   shellQuote: False
 action2:
   type: string
   default: "sort"
   inputBinding:
    position: 6
    shellQuote: False
 threads2:
   type: int
   inputBinding:
     prefix: -@
     position: 7
     shellQuote: False
 outfilename:
   type: string
   inputBinding:
     prefix: -o
     position: 8
     shellQuote: False

outputs:
  example_out:
   type: File
   outputBinding:
    glob: $(inputs.outfilename)
