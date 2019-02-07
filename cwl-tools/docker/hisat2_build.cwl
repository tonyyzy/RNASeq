#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: hisat2-build
hints:
   DockerRequirement:
      dockerPull: quay.io/biocontainers/hisat2:2.1.0--py27h2d50403_2 
stdout: log.txt
inputs:
   reference:
      type: File
      inputBinding:
         position: 1
         prefix: -f
   basename:
      type: string
      inputBinding:
         position: 2
   threads:
      type: int
      inputBinding:
         prefix: -p

outputs:
   ht:
      type: File[]
      outputBinding:
         glob: "*"
   log:
      type: stdout
