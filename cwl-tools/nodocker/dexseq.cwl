r/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript

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

outputs:
   example_out:
      type: File
      outputBinding:
         glob: "*"
