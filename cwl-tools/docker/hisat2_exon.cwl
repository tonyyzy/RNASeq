#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: hisat2_extract_exons.py
requirements:
   InlineJavascriptRequirement: {}
hints:
   DockerRequirement:
      dockerPull: quay.io/biocontainers/hisat2:2.1.0--py27h2d50403_2 
stdout: $(inputs.gtf.basename + ".exons")
inputs:
   gtf:
      type: File
      inputBinding:
         position: 1

outputs:
   exons:
      type: stdout