#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand:
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
arguments:
  - position: 0
    shellQuote: False
    valueFrom: $("mkdir " + inputs.outFileNamePrefix + "_STARAligner" + " && cd " + inputs.outFileNamePrefix + "_STARAligner" + " && STAR ")
  - position: 5
    valueFrom: |
      ${
        if (inputs.readFilesIn[0].nameext == ".gz"){
          return "--readFilesCommand gunzip -c";}
          return "";
      }

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/star:2.6.0c--0

inputs:
  Threads:
    type: int
    inputBinding:
      position: 1
      prefix: --runThreadN
  genomeDir:
    type: Directory
    inputBinding:
      position: 2
      prefix: --genomeDir
  readFilesIn:
    type: File[]
    inputBinding:
      position: 3
      prefix: --readFilesIn
  outFileNamePrefix:
    type: string
    inputBinding:
      position: 4
      prefix: --outFileNamePrefix
  XSTag:
    type: string?
    default: intronMotif
    inputBinding:
      position: 5
      prefix: --outSAMstrandField

outputs:
  star_read_out:
    type: Directory
    outputBinding:
      glob: $(inputs.outFileNamePrefix + "_STARAligner")

  sam_output:
    type: File
    outputBinding:
      glob: $(inputs.outFileNamePrefix + "_STARAligner/*.sam")
