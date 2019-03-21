#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript

requirements:
   DockerRequirement:
      dockerPull: filipejesus/dexseq:3.8

inputs:
   input_script:
      type: File
      inputBinding:
         position: 1
   counts_matrix:
      type: Directory
      inputBinding:
         position: 2
         prefix: --count_matrix_dir
   metadata:
      type: File
      inputBinding:
         position: 4
         prefix: --metadata
   threads:
      type: int?
      inputBinding:
        position: 5
        prefix: --threads

outputs:
   dexseq_out:
      type: File[]
      outputBinding:
         glob: "*.csv"
