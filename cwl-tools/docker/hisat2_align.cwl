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
   -  position: 0
      shellQuote: False
      valueFrom: $("mkdir " + inputs.sam_name.split(".")[0] + " && cd " + inputs.sam_name.split(".")[0] + " && hisat2")
inputs:
   input_type:
      type: string
      default: "-q"
      inputBinding:
         position: 1
   index_directory:
      type: Directory
      inputBinding:
         position: 2
         prefix: "-x"
         valueFrom: $(inputs.index_directory.path + "/" +
            inputs.index_directory.listing[0].nameroot.split(".").slice(0,-1).join("."))
   first_pair:
      type: File?
      inputBinding:
         position: 3
         prefix: "-1"
   second_pair:
      type: File?
      inputBinding:
         position: 4
         prefix: "-2"
   single_file:
      type: File[]?
      inputBinding:
         position: 3
         prefix: "-U"
   sra_acc:
      type: string?
      inputBinding:
         position: 3
         prefix: "--sra-acc"
   sam_name:
      type: string
      inputBinding: 
         prefix: -S
   threads:
      type: int
      inputBinding:
         prefix: -p
   XSTag:
      type: string?
      inputBinding:
         position: 5
   log:
      type: string
      default: "log.txt"
      inputBinding:
         prefix: --summary-file

outputs:
   hisat2_align_out:
      type: Directory
      outputBinding:
         glob: $(inputs.sam_name.split(".")[0])
   
   sam_output:
      type: File
      outputBinding:
         glob: $(inputs.sam_name.split(".")[0] + "/" + inputs.sam_name)

