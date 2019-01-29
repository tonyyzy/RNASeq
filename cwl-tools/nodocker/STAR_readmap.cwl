#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand:
requirements:
  - class: ShellCommandRequirement
arguments: ["mkdir", $(inputs.outFileNamePrefix), "&&", "cd", $(inputs.outFileNamePrefix), "&&", "STAR"]

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
  out:
    type: Directory
    outputBinding:
      glob: $(inputs.outFileNamePrefix)
