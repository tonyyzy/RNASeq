#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: hisat2_extract_splice_sites.py
hints:
   DockerRequirement:
      dockerPull: quay.io/biocontainers/hisat2:2.1.0--py27h2d50403_2 
stdout: $(inputs.gtf.basename + ".splice_sites")
inputs:
   gtf:
      type: File
      inputBinding:
         position: 1

outputs:
   splice_sites:
      type: stdout