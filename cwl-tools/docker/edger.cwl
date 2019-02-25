#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement

hints:
    DockerRequirement:
        dockerPull: quay.io/biocontainers/bioconductor-edger:3.24.1--r351hf484d3e_0

inputs:
    script:
        type: File
        inputBinding:
            position: 0
    condition:
        type: string
        inputBinding:
            position: 1
            prefix: --condition
    counts:
        type: File
        inputBinding:
            position: 2
            prefix: --counts
    metadata:
        type: File
        inputBinding:
            position: 3
            prefix: --metadata

outputs:
    output:
      type: File
      outputBinding:
        glob: $("DGE_reds.csv")
