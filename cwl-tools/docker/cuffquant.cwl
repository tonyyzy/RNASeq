#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: cuffquant

requirements:
    InlineJavascriptRequirement: {}

hints:
    DockerRequirement:
        dockerPull: machalen/cufflinksdocker:latest

inputs:
    output:
        type: string
        inputBinding:
            position: 1
            prefix: -o
    threads:
        type: int
        inputBinding:
            position: 2
            prefix: -p
    bias_correction:
        type: File?
        inputBinding:
            position: 3
            prefix: -b
    multi_read_correct:
        type: string?
        inputBinding:
            position: 4
            prefix: -u
    libType:
        type: string?
        inputBinding:
            prefix: –library-type
            position: 5
    frag_mean_len:
        type: string?
        inputBinding:
            prefix: -m
            position: 6
    frad_sd_len:
        type: string?
        inputBinding:
            prefix: -string
            position: 7
    remove_len_correction:
        type: string?
        inputBinding:
            prefix: –no-length-correction
            position: 8
    merged_gtf:
        type: File
        inputBinding:
            position: 9
    alignment_file:
        type: File
        inputBinding:
            position: 10

outputs:
   output:
      type: Directory
      outputBinding:
         glob: $(inputs.output)
   cxb_out:
      type: File
      outputBinding:
         glob: $(inputs.output+"/abundances.cxb")
         outputEval: |
          ${
            self[0].basename = inputs.output + '.cxb';
            return self[0]
          }
