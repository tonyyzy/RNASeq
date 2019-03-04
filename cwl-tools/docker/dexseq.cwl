#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript
hints:
   DockerRequirement:
      dockerPull: machalen/dexseq

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
   gff:
      type: Directory
      inputBinding:
         position: 3
         prefix: --gff_file_dir
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
   output:
      type: File[]
      outputBinding:
         glob: "*.csv"
