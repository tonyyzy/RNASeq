#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: python2
<<<<<<< HEAD
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
=======
inputs:
 program:
  type: string
  inputBinding:
    position: 1
 gtfDir:
  type: Directory[]
  inputBinding:
>>>>>>> 41d1100c4f6a457b46f750e5d2f93194d78e0034
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
