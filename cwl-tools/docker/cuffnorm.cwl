#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: mkdir

requirements:
    ShellCommandRequirement: {}

hints:
    DockerRequirement:
        dockerPull: machalen/cufflinksdocker:latest

arguments:
  - position: -4
    valueFrom: $(inputs.output)
  - position: -3
    prefix: "&&"
    valueFrom: "cuffnorm"
    shellQuote: false

inputs:
  output:
     type: string
     inputBinding:
        position: 1
        prefix: -o
  threads:
     type: int
     inputBinding:
        position: 2
        prefix: -p
  merged_gtf:
     type: File
     inputBinding:
        position: 3
  condition1_files:
     type: File[]
     inputBinding:
        itemSeparator: ","
        separate: false
        position: 4
        prefix: ""
  condition2_files:
     type: File[]
     inputBinding:
        itemSeparator: ","
        separate: false
        position: 5
        prefix: ""

outputs:
   cuffnorm_out:
    type: Directory
    outputBinding:
      glob: $(inputs.output)
