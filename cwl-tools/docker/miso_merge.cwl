cwlVersion: v1.0
class: CommandLineTool
baseCommand:

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  DockerRequirement:
     dockerPull: genomicpariscentre/miso:0.5.3

arguments:
   - position: 1
     valueFrom: $("samtools merge " + inputs.output)
     shellQuote: False

inputs:
   bam:
      type: File[]
      inputBinding:
         position: 2
   output:
      type: string

outputs:
   miso_out:
      type: File
      outputBinding:
         glob: $(inputs.output)
