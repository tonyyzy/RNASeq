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
   sam:
      type: File?
      inputBinding:
         position: 3
   gff:
      type: File?
      inputBinding:
         position: 2
   gff_name:
      type: string?
      inputBinding:
         position: 3
   outname:
      type: string?
      inputBinding:
         position: 4


outputs:
   output:
      type: File
      outputBinding:
         glob: "*"
