#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand:
hints:
   DockerRequirement:
      dockerPull: filipejesus/miso:latest
requirements:
   InlineJavascriptRequirement: {}
   ShellCommandRequirement: {}

arguments:
   - position: 0
     valueFrom: cp $(inputs.bam.path) $(runtime.outdir)
     shellQuote: False
   - position: 1
     valueFrom: "&&"
   - position: 2
     valueFrom: cp -r $(inputs.index.path) $(runtime.outdir)
     shellQuote: False
   - position: 3
     valueFrom: "&&"
   - position: 4
     valueFrom: sh /complete_run.sh
     shellQuote: False
   - position: 7
     valueFrom: $(runtime.outdir+"/"+inputs.index.basename)
     shellQuote: False
   - position: 8
     valueFrom: $(runtime.outdir+"/"+inputs.bam.basename)
     shellQuote: False
   - position: 10
     valueFrom: $(runtime.outdir+"/"+inputs.out_dir)
     shellQuote: False


inputs:
   cores:
      type: string
      inputBinding:
         position: 5
   lib_type:
      type: string
      inputBinding:
         position: 6
   index:
      type: Directory
   bam:
      type: File
   read_len:
      type: string
      inputBinding:
         position: 9
   out_dir:
      type: string
   annotation_file:
      type: File?
      inputBinding:
         position: 11
   exon_size:
      type: string?
      inputBinding:
         position: 12

outputs:
   output:
      type: Directory
      outputBinding:
         glob: $(inputs.out_dir)
