#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: index_gff
hints:
   DockerRequirement:
      dockerPull: filipejesus/miso 

inputs:
   annotation:
      type: File
      inputBinding:
         position: 1
         prefix: --index
   index_dir:
      type: string
      inputBinding:
         position: 2

outputs:
   output:
      type: Directory
      outputBinding:
         glob: $(inputs.index_dir)
