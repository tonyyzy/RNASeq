#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: 
hints:
  ShellCommandRequirement: {}
  DockerRequirement:
    dockerPull: quay.io/biocontainers/hisat2:2.1.0--py27h2d50403_2
requirements:
  InlineJavascriptRequirement: {}
arguments:
  - position: 1
    valueFrom: "mkdir"
  - position: 2
    valueFrom: $(inputs.sam_name.split('.')[0])
  - position: 3
    shellQuote: False
    valueFrom: '&& cd'
  - position: 4
    valueFrom: $(inputs.sam_name.split('.')[0])
  - position: 5
    shellQuote: False
    valueFrom: '&& hisat2'

inputs:
  input_type:
    type: string
    default: "-q"
    inputBinding:
      position: 6
  index_directory:
    type: Directory
    inputBinding:
      position: 7
      prefix: "-x"
      valueFrom: "${return inputs.index_directory.path + '/' 
                  + inputs.index_directory.listing[0].nameroot.split('.').slice(0,-1).join('.')}"
  first_pair:
    type: File?
    inputBinding:
      position: 8
      prefix: "-1"
  second_pair:
    type: File?
    inputBinding:
      position: 9
      prefix: "-2"
  single_file:
    type: File[]?
    inputBinding:
      position: 10
      prefix: -U
  sra_acc:
    type: string?
    inputBinding:
      position: 11
      prefix: --sra-acc
  sam_name:
    type: string
    inputBinding:
      position: 12
      prefix: -S
  threads:
    type: int
    inputBinding:
      position: 13
      prefix: -p
  XSTag:
    type: string?
    default: --dta-cufflinks
    inputBinding:
      position: 14
  log:
    type: string
    default: "log.txt"
    inputBinding:
      position: 15
      prefix: --summary-file

outputs:
  hisat2_align_out:
    type: Directory
    outputBinding:
      glob: $(inputs.sam_name.split('.')[0])
  
  sam_output:
    type: File
    outputBinding:
      glob: $(inputs.sam_name.split('.')[0] + '/' + inputs.sam_name)

