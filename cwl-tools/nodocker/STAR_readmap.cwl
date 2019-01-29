#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand:
requirements:
  - class: ShellCommandRequirement
arguments: ["mkdir", $(inputs.outFileNamePrefix), "&&", "cd", $(inputs.outFileNamePrefix), "&&", "STAR"]

inputs:
  Threads:
    type: string
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
#  sam_output:
#    type: File
#    outputBinding:
#      glob: "*.sam"
#  log_outputs:
#    type: File[]
#    outputBinding:
#      glob: "*.out"
#  tab_output:
#    type: File
#    outputBinding:
#      glob: "*.tab"
