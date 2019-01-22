#!user/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: dexseq_prepare_annotation.py 

inputs:
   GTF_input:
      type: File
      inputBinding:
         position: 1
   GFF_output_name:
      type: string
      inputBinding:
         position: 2

outputs:
   outout:
      type: File[]
      outputBinding:
         glob: "*"
