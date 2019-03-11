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
     valueFrom: $("samtools merge " + inputs.output[0])
     shellQuote: False
   - prefix: "&&"
     position: 3
     valueFrom: $("samtools merge " + inputs.output[1])
     shellQuote: False

inputs:
   condition1bam:
      type: File[]
      inputBinding:
         position: 2
   condition2bam:
      type: File[]
      inputBinding:
         position: 4
   output:
      type: string[]

outputs:
   condition1:
      type: File
      outputBinding:
         glob: $(inputs.output[0])
   condition2:
      type: File
      outputBinding:
         glob: $(inputs.output[1])
