#!usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand:
hints:
   DockerRequirement:
      dockerPull: filipejesus/miso
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
     valueFrom: cp $(inputs.bam_bai.path) $(runtime.outdir)
     shellQuote: False
   - position: 3
     valueFrom: "&&"
   - position: 4
     valueFrom: cp -r $(inputs.index.path) $(runtime.outdir)
     shellQuote: False
   - position: 5
     valueFrom: "&&"
   - position: 6
     valueFrom: sh /complete_run.sh
     shellQuote: False
   - position: 9
     valueFrom: $(runtime.outdir+"/"+inputs.index.basename)
     shellQuote: False
   - position: 10
     valueFrom: $(runtime.outdir+"/"+inputs.bam.basename)
     shellQuote: False
   - position: 12
     valueFrom: $(runtime.outdir+"/"+inputs.out_dir)
     shellQuote: False


inputs:
   cores:
      type: string
      inputBinding:
         position: 7
   lib_type:
      type: string
      inputBinding:
         position: 8
   index:
      type: Directory
   bam:
      type: File
   bam_bai:
      type: File
   read_len:
      type: string
      inputBinding:
         position: 11
   out_dir:
      type: string
   annotation_file:
      type: File?
      inputBinding:
         position: 13
   exon_size:
      type: string?
      inputBinding:
         position: 14

outputs:
   output:
      type: Directory
      outputBinding:
         glob: $(inputs.out_dir)
