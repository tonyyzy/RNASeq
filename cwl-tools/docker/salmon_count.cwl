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
   salmon_dir:
      type: Directory
      inputBinding:
         position: 4
         prefix: --salmon_dir

outputs:
   output:
      type: File
      outputBinding:
         glob: "*gene_count_matrix*"
