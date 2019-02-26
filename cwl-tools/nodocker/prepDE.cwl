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
        for (var i = 0; i < inputs.stringtie_out.length; i ++) {
          paths.push({
          "class": "Directory",
          "basename": inputs.stringtie_out[i].nameroot,
          "listing": [inputs.stringtie_out[i]]
          });
        }return paths;
      }
arguments:
  - prefix: -i
    position: 2
    valueFrom: $(runtime.outdir)
inputs:
  input_script:
    type: File
    inputBinding:
      position: 1
  stringtie_out:
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
