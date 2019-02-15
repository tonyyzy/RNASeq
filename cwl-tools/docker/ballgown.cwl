#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: Rscript

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bioconductor-ballgown:2.14.0--r351_0

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: | 
      ${
        var paths = [];
        for (var i = 0; i < inputs.tablemaker_output.length; i ++) {
          paths.push({
            "class": "Directory",
            "basename": inputs.tablemaker_output[i].basename,
            "listing": inputs.tablemaker_output[i].listing
          });
        }return paths;
      }

arguments:
  - prefix: --data_dir
    position: 1
    valueFrom: $(runtime.outdir)

inputs:
  script:
    type: File
    inputBinding:
      position: 0
  tablemaker_output:
    type: Directory[]
  metadata:
    type: File
    inputBinding:
      position: 2
      prefix: --metadata
  condition:
    type: string
    inputBinding:
      position: 3
      prefix: --condition

outputs:
  gene_matrix:
    type: File
    outputBinding:
      glob: $("gene_count_matrix.csv")
  transcript_matrix:
    type: File
    outputBinding:
      glob: $("transcript_count_matrix.csv")
