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
   index_name:
      type: string
      inputBinding:
         position: 2
         prefix: -i
   index_type:
      type: string
      inputBinding:
         position: 3
         prefix: --type
   cores:
      type: string
      inputBinding:
         position: 4
         prefix: -p

outputs:
   output:
      type: Directory
      outputBinding:
         glob: $(inputs.index_name)
