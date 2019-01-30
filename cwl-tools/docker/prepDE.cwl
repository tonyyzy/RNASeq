#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: python2
hints:
   DockerRequirement:
      dockerPull: python:2.7.15-slim

requirements:
  ShellCommandRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: '{"inputs": $(inputs.gtfs)}'
        entryname: inputs.json

arguments:
  - prefix: -c
    position: 2
    valueFrom: |
      import json
      with open("inputs.json") as file:
          inputs = json.load(file)
      with open("sample_lst.txt", "w") as txt:
          for i in range(len(inputs["inputs"])):
              txt.write(inputs["inputs"][i]["nameroot"]
                        + " "
                        + inputs["inputs"][i]["path"]
                        + "\n")
  - prefix: "&&"
    position: 3
    valueFrom: "python2"
  - prefix: -i
    position: 5
    valueFrom: sample_lst.txt
inputs:
  program:
    type: File
    inputBinding:
      position: 4
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
