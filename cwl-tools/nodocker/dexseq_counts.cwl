#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: dexseq_count.py

inputs:
   gff_input:
      type: File
      inputBinding:
         position: 1
   sam_input:
      type: File
      inputBinding:
         position: 2
   count_file_name:
      type: string
      inputBinding:
         position: 3

outputs:
   output:
      type: File[]
      outputBinding:
         glob: "*"
