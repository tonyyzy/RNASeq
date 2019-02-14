#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: cufflinks

hints:
   DockerRequirement:
      dockerPull: machalen/cufflinksdocker:latest

inputs:
   gtf:
      type: File
      inputBinding:
         position: 1
         prefix: -G
   output_dir:
      type: string
      inputBinding:
         position: 2
         prefix: -o
   threads:
      type: string
      inputBinding:
         position: 3
         prefix: -p
   BiasCorrecting:
      type: File?
      inputBinding:
         position: 4
         prefix: -b
   libType:
      type: string?
      inputBinding:
         position: 5
         prefix: –library-type
   libNorm:
      type: string?
      inputBinding:
         position: 6
         prefix: –library-norm-method
   multi_read_correction:
      type: string?
      inputBinding:
         position: 7
         prefix: -u
   alignment_file:
      type: File
      inputBinding:
         position: 8

outputs:
   output:
      type: Directory
      outputBinding:
         glob: $(inputs.output_dir)