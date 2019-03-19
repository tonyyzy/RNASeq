#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand:

requirements:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/stringtie:1.3.0--hd28b015_2
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}

arguments:
 - position: -4
   valueFrom: "stringtie"
   shellQuote: False
 - position: 4
   prefix: -o
   valueFrom: $(inputs.output + "/" + inputs.output + ".gtf")

inputs:
 bam:
  type: File
  inputBinding:
   position: 1
   prefix: -eB
 threads:
  type: int
  inputBinding:
   position: 2
   prefix: -p
 gtf:
  type: File
  inputBinding:
   position: 3
   prefix: -G
 output:
  type: string

outputs:
 stringtie_out:
  type: Directory
  outputBinding:
   glob: $(inputs.output)
