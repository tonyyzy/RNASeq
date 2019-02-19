cwlVersion: v1.0
class: CommandLineTool
baseCommand: salmon
hints:
   DockerRequirement:
      dockerPull: combinelab/salmon:0.12.0
requirements:
   InlineJavascriptRequirement: {}

arguments:
   - position: -1
     valueFrom: "quant"
   - position: 0
     valueFrom: "--validateMappings"
   - position: 4
     valueFrom: "A"
     prefix: "-l"
inputs:
   index:
      type: Directory
      inputBinding:
         position: 1
         prefix: "-i"
   out_dir:
      type: string
      inputBinding:
         position: 2
         prefix: "-o"
   cores:
      type: int
      inputBinding:
         position: 3
         prefix: "-p"
   first_end_fastq:
      type: File?
      inputBinding:
         position: 5
         prefix: "-1"
   second_end_fastq:
      type: File?
      inputBinding:
         position: 6
         prefix: "-2"
   single_fastq:
      type: File[]?
      inputBinding:
         position: 5
         prefix: "-r"

outputs:
   output:
      type: Directory
      outputBinding:
         glob: $(inputs.out_dir)
