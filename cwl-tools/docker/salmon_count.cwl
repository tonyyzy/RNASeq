#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript
hints:
   DockerRequirement:
      dockerPull: filipejesus/tximport

inputs:
   input_script:
      type: File
      inputBinding:
         position: 1
   gtf:
      type: File
      inputBinding:
         position: 2
         prefix: --gtf
   metadata:
      type: File
      inputBinding:
         position: 3
         prefix: --metadata
   quant_results:
      type: Directory
      inputBinding:
         position: 4
         prefix: --salmon_dir

outputs:
   count:
      type: File
      outputBinding:
         glob: "*gene_count_matrix*"
   length:
      type: File
      outputBinding:
         glob: "*gene_length_matrix*"
   abundance:
      type: File
      outputBinding:
         glob: "*gene_abundance_matrix*"
