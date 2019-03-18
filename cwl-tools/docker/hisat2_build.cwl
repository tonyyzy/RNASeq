#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: hisat2-build
hints:
   DockerRequirement:
      dockerPull: quay.io/biocontainers/hisat2:2.1.0--py27h2d50403_2 
stdout: hisat2_build.log
stderr: hisat2_build.log
inputs:
   threads:
      type: int
      inputBinding:
         position: 1
         prefix: -p
   splice_sites:
      type: File?
      inputBinding:
         position: 2
         prefix: --ss
   exon:
      type: File?
      inputBinding:
         position: 3
         prefix: --exon
   fasta:
      type: File[]
      inputBinding:
         itemSeparator: ","
         separate: True
         position: 4
         prefix: -f
         shellQuote: False
   output:
      type: string
      inputBinding:
         position: 5

outputs:
   ht:
      type: File[]
      outputBinding:
         glob: "*.ht2*"
   log:
      type: stdout
