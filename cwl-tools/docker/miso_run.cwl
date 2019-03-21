#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: sh

hints:
   DockerRequirement:
      dockerPull: genomicpariscentre/miso:0.5.3

inputs:
   input_script:
      type: File
      inputBinding:
         position: 1
   gff:
      type: File
      inputBinding:
         position: 6
   bam:
      type: File
      inputBinding:
         position: 2
   lib_type:
      type: string
      inputBinding:
         position: 4
   min_exon_size:
      type: string
      default: "50"
      inputBinding:
         position: 3
   output:
      type: string
      inputBinding:
         position: 5

outputs:
   miso_out:
      type: Directory
      outputBinding:
         glob: $(inputs.output)
