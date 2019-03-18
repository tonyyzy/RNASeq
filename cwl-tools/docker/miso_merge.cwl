cwlVersion: v1.0
class: CommandLineTool
baseCommand:
hints:
   DockerRequirement:
      dockerPull: genomicpariscentre/miso

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}

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
