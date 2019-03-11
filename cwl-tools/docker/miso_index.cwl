#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: index_gff
hints:
   DockerRequirement:
      dockerPull: filipejesus/miso

inputs:
   gtf:
      type: File
      inputBinding:
         position: 1
         prefix: --index
   output:
      type: string
      inputBinding:
         position: 2

outputs:
   miso_out:
      type: Directory
      outputBinding:
         glob: $(inputs.output)
