#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: perl

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  DockerRequirement:
     dockerPull: genomicpariscentre/miso:0.5.3

arguments:
   - position: 3
     prefix: "|"
     valueFrom: $("tee " + runtime.outdir + "/" + inputs.output)
     shellQuote: False

inputs:
   perl_input:
      type: File
      inputBinding:
         position: 1
   gtf:
      type: File
      inputBinding:
         position: 2
   output:
      type: string

outputs:
   miso_out:
      type: File
      outputBinding:
         glob: $(inputs.output)
