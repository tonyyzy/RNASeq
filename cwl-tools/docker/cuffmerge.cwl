#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: cuffmerge

hints:
    DockerRequirement:
        dockerPull: machalen/cufflinksdocker:latest

requirements:
    ShellCommandRequirement: {}
    InitialWorkDirRequirement:
        listing:
            - entry: '{"inputs": $(inputs.cufflinks_output)'
              entryname: inputs.json

arguments:
    - position: 4
      valueFrom: |
        import json
        with open("inputs.json") as file:
            inputs = json.load(file)
        with open("assembly_GTF_list.txt") as txt:
            for i in range(len(inputs["inputs"])):
                txt.write(inputs["inputs"][i]["path"] + " \n")

inputs:
    output:
        type: string
        inputBinding:
            position: 1
            prefix: -o
    gtf:
        type: File
        inputBinding:
            position: 2
            prefix: -g
    threads:
        type: File
        inputBinding:
            position: 3
            prefix: -p
    cufflinks_output:
        type: File[]

outputs:
    output:
        type: File
        outputBinding:
            glob: $(inputs.output)