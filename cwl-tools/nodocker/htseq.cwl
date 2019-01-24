r/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: python2.7

inputs:
   input_script:
      type: File
      inputBinding:
         position: 1 
   gtf:
      type: File?
      inputBinding:
         position: 2
   sam:
      type: File?
      inputBinding:
         position: 3
   gff:
      type: File?
      inputBinding:
         position: 2
   gff_name:
      type: string?
      inputBinding:
         position: 3
   outname:
      type: string?
      inputBinding:
         position: 4


outputs:
   example_out:
      type: File
      outputBinding:
         glob: "*"
