#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript
hints:
   DockerRequirement:
      dockerPull: filipejesus/fgsea

inputs:
   input_scripts:
      type: File
      inputBinding:
         position: 1
   de_results:
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
   example_out:
      type: File
      outputBinding:
         glob: "gsea_res.csv"
