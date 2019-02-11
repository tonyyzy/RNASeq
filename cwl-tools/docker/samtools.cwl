#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: genomicpariscentre/samtools
baseCommand: samtools

inputs:
  action:
    type: string
    default: "sort"
    inputBinding:
      position: 1
  sortby:
    type: string
    default:
    inputBinding:
      position: 2
  threads:
    type: int
    inputBinding:
      prefix: -@
      position: 3
  outfilename:
    type: string
    inputBinding:
      prefix: -o
      position: 4
  samfile:
    type: File
    inputBinding:
      position: 5

outputs:
  samtools_out:
    type: File
    outputBinding:
      glob: $(inputs.outfilename)
