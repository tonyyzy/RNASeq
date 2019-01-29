#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: python2


requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: |
      ${
        var paths = [];
        for (var i = 0; i < inputs.gtfs.length; i ++) {
          paths.push({
          "class": "Directory",
          "basename": inputs.gtfs[i].nameroot,
          "listing": [inputs.gtfs[i]]
          });
        }return paths;
      }
arguments:
  - prefix: -i
    position: 2
    valueFrom: $(runtime.outdir)
inputs:
  program:
    type: File
    inputBinding:
      position: 1
  gtfs:
    type: File[]

outputs:
  gene_output:
    type: File
    outputBinding:
      glob: "*gene_count_matrix*"
  transcript_output:
    type: File
    outputBinding:
      glob: "*transcript_count_matrix*"
