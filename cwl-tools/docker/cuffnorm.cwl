#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: cuffnorm

arguments:
  - position: -4
    valueFrom: $(inputs.output)
  - position: -3
    prefix: "&&"
    valueFrom: "cuffdiff"
    shellQuote: false

inputs:
  output:
     type: string
     inputBinding:
        position: 1
        prefix: -o
  threads:
     type: string
     inputBinding:
        position: 2
        prefix: -p`
  gtf_file:
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
        arguments:
           - position: -4
             valueFrom: $(inputs.output)
           - position: -3
             prefix: "&&"
             valueFrom: "cuffdiff"
             shellQuote: false
  condition2_files:
     type: File[]
     inputBinding:
        itemSeparator: ","
        separate: false
        position: 5
        prefix: ""

outputs:
   output:
    type: Directory
    outputBinding:
      glob: $(inputs.output)
