cwlVersion: v1.0
class: CommandLineTool
baseCommand: salmon
hints:
   DockerRequirement:
      dockerPull: combinelab/salmon:0.12.0
requirements:
   InlineJavascriptRequirement: {}

arguments:
   - position: 0
     valueFrom: "index"
inputs:
   fasta:
      type: File
      inputBinding:
         position: 1
         prefix: -t
   output:
      type: string
      inputBinding:
         position: 2
         prefix: -i
   index_type:
      type: string
      inputBinding:
         position: 3
         prefix: --type
   threads:
      type: int
      inputBinding:
         position: 4
         prefix: -p

outputs:
   salmon_out:
      type: Directory
      outputBinding:
         glob: $(inputs.output)
