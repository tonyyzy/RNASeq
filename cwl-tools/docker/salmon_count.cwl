#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript

requirements:
   DockerRequirement:
      dockerPull: filipejesus/tximport:3.8

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
   gene_count_output:
      type: File
      outputBinding:
         glob: "*gene_count_matrix*"
   gene_length_output:
      type: File
      outputBinding:
         glob: "*gene_length_matrix*"
   gene_abundance_output:
      type: File
      outputBinding:
         glob: "*gene_abundance_matrix*"
