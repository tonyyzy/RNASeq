cwlVersion: v1.0
class: CommandLineTool
baseCommand: python2

hints:
    DockerRequirement:
        dockerPull: filipejesus/cufflinks:latest

requirements:
    ShellCommandRequirement: {}
    InlineJavascriptRequirement: {}
    InitialWorkDirRequirement:
        listing:
            - entry: '{"inputs": $(inputs.cufflinks_output)}'
              entryname: inputs.json

arguments:
    - position: -3
      prefix: -c
      valueFrom: |
        import json
        with open("inputs.json") as file:
            inputs = json.load(file)
        with open("assembly_GTF_list.txt", "w") as txt:
            for i in range(len(inputs["inputs"])):
                txt.write(inputs["inputs"][i]["path"]
                    + "\n")
    - prefix: "&&"
      position: -2
      shellQuote: False
      valueFrom: $("cp "+ inputs.fasta.path + " " + runtime.outdir)
    - prefix: "&&"
      position: -1
      valueFrom: "cuffmerge"
    - position: 5
      valueFrom: assembly_GTF_list.txt
    - position: 4
      prefix: -s
      valueFrom: $(runtime.outdir+"/"+inputs.fasta.basename)
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
        type: int
        inputBinding:
            position: 3
            prefix: -p
    fasta:
        type: File
#        inputBinding:
#            position: 4
#            prefix: -s
    cufflinks_output:
        type: File[]

outputs:
    cuffmerge_out:
        type: Directory
        outputBinding:
            glob: $(inputs.output)
    merged_gtf:
        type: File
        outputBinding:
            glob: $(inputs.output+"/merged.gtf")
