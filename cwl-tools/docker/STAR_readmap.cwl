#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand:
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
arguments: ["mkdir", $(inputs.outFileNamePrefix + "_STARAligner"), "&&", "cd", $(inputs.outFileNamePrefix + "_STARAligner"), "&&", "STAR"]

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/star:2.6.0c--0

inputs:
  Threads:
    type: int
    inputBinding:
      prefix: --runThreadN
  genomeDir:
    type: Directory
    inputBinding:
      prefix: --genomeDir
  readFilesIn:
    type: File[]
    inputBinding:
      prefix: --readFilesIn
  outFileNamePrefix:
    type: string?
    inputBinding:
      prefix: --outFileNamePrefix

outputs:
  star_read_out:
    type: Directory
    outputBinding:
      glob: $(inputs.outFileNamePrefix + "_STARAligner")

  sam_output:
    type: File
    outputBinding:
      glob: $(inputs.outFileNamePrefix + "_STARAligner/*.sam")
