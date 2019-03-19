#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: compare_miso

requirements:
   DockerRequirement:
      dockerPull: genomicpariscentre/miso

inputs:
   group1:
      type: Directory
      inputBinding:
         position: 1
         prefix: "--compare-samples"
   group2:
      type: Directory
      inputBinding:
         position: 2
   output:
      type: string
      inputBinding:
         position: 3

outputs:
   miso_out:
      type: Directory
      outputBinding:
         glob: $(inputs.output)
