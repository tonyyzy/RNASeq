#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: mkdir

requirements:
    ShellCommandRequirement: {}
    InlineJavascriptRequirement: {}
    DockerRequirement:
        dockerPull: filipejesus/cufflinks:2.2.1

arguments:
  - position: -4
    valueFrom: $(inputs.output)
  - position: -3
    prefix: "&&"
    valueFrom: "cuffnorm"
    shellQuote: False
  - position: 6
    prefix: "&&"
    valueFrom: $("python")
    shellQuote: False
  - position: 8
    valueFrom: $(inputs.output + "/genes.count_table")
  - position: 9
    valueFrom: $(inputs.output + "/genes.attr_table")


inputs:
  input_script:
     type: File
     inputBinding:
        position: 7
  output:
     type: string
     inputBinding:
        position: 1
        prefix: -o
  metadata:
     type: File
     inputBinding:
        position: 7
  threads:
     type: int
     inputBinding:
        position: 2
        prefix: -p
  label:
     type: string[]
     inputBinding:
        position: 10
        itemSeparator: ","
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
