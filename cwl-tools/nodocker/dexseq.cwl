#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: DEXSeq.R

inputs:
   matrix:
      type: string
      inputBinding:
         prefix: --count_matrix_dir
   gff_file:
      type: File
      inputBinding:
         prefix: --gff_file_dir
   metadata:
      type: File
      inputBinding:
         prefix: --metadata

outputs:
   output:
      type: File[]
      outputBinding:
      glob: "*"
