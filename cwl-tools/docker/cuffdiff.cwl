#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: cuffdiff

hints:
   DockerRequirement:
       dockerPull: machalen/cufflinksdocker:latest

inputs:
   output:
      type: string
      inputBinding:
         position: 1
         prefix: -o
   threads:
      type: string
      inputBinding:
         position: 2
         prefix: -p
   label:
      type: string[]
      inputBinding:
         position: 3
         itemSeparator: ","
         prefix: -L
   use_all:
      type: string?
      inputBinding:
         position: 4
         prefix: –total-hits-norm
   compatible_only:
      type: string?
      inputBinding:
         position: 5
         prefix: –compatible-hits-norm
   bias_correction:
      type: File?
      inputBinding:
         position: 6
         prefix: -b
   multi_read_correct:
      type: string?
      inputBinding:
         position: 7
         prefix: -u
   FDR:
      type: string
      inputBinding:
         position: 8
         prefix: –FDR
   libType:
      type: string
      inputBinding:
         position: 9
         prefix: –library-type
   libNorm:
      type: string
      inputBinding:
         position: 10
         prefix: –library-norm-method
   gtf_file:
      type: File
      inputBinding:
         position: 11
   condition1_files:
      type: File[]
      inputBinding:
         itemSeparator: ","
         separate: false
         position: 12
   condition2_files:
      type: File[]
      inputBinding:
         itemSeparator: ","
         separate: false
         position: 13

outputs:
   output:
      type: Directory
      outputBinding:
         glob: $(inputs.output)
