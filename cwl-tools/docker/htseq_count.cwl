#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: python
hints:
   DockerRequirement:
      dockerPull: genomicpariscentre/htseq:0.11.0

inputs:
   input_script:
      type: File
      inputBinding:
         position: 1
   pairedend:
      type: string
      inputBinding:
         position: 2
         prefix: -p
   stranded:
      type: string
      inputBinding:
         position: 3
         prefix: -s
   input_format:
      type: string
      inputBinding:
         position: 4
         prefix: -f
   sorted_by:
      type: string
      inputBinding:
         position: 5
         prefix: -r
   gff:
      type: File?
      inputBinding:
         position: 6
   sam:
      type: File?
      inputBinding:
         position: 7
   outname:
      type: string?
      inputBinding:
         position: 8


outputs:
   output:
      type: File
      outputBinding:
         glob: "*"
