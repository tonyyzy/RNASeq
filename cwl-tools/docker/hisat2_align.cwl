#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: hisat2
hints:
   DockerRequirement:
      dockerPull: quay.io/biocontainers/hisat2:2.1.0--py27h2d50403_2
requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.index_directory)

inputs:
   input_type:
      type: string
      default: "-q"
      inputBinding:
         position: 1
   index_directory:
      type: Directory
   index_basename:
      type: string
      inputBinding:
         position: 2
         prefix: "-x"
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
      type: File?
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

outputs:
   output:
      type: File[]
      outputBinding:
         glob: "*.sam"
