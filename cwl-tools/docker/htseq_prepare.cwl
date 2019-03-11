#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: python
hints:
   DockerRequirement:
      dockerPull: genomicpariscentre/htseq

inputs:
   input_script:
      type: File
      inputBinding:
         position: 1
   gtf:
      type: File?
      inputBinding:
         position: 2
   gff_name:
      type: string?
      inputBinding:
         position: 3


outputs:
   ht_prep_out:
      type: File
      outputBinding:
         glob: "*"
