#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript
hints:
   DockerRequirement:
      dockerPull: filipejesus/fgsea

inputs:
   input_script:
      type: File
      inputBinding:
         position: 1
   de_res:
      type: File
      inputBinding:
         position: 2
         prefix: --de_res
   gene_set:
      type: File
      inputBinding:
         position: 3
         prefix: --gene_set

outputs:
   gsea_out:
      type: File
      outputBinding:
         glob: "gsea_res.csv"
