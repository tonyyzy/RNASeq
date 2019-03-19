cwlVersion: v1.0
class: CommandLineTool
baseCommand: salmon

requirements:
   DockerRequirement:
      dockerPull: combinelab/salmon:0.12.0
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
   index_directory:
      type: Directory
      inputBinding:
         position: 1
         prefix: "-i"
   output:
      type: string
      inputBinding:
         position: 2
         prefix: "-o"
   threads:
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
   salmon_out:
      type: Directory
      outputBinding:
         glob: $(inputs.output)
