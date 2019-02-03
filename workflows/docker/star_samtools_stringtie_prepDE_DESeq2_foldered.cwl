#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
requirements:
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}

inputs:
  Threads: int
  genomeDir: Directory
  readFilesIn_1: File[]
  readFilesIn_2: File[]
  outFileNamePrefix_1: string
  outFileNamePrefix_2: string
  outfilename_samtools_1: string
  outfilename_samtools_2: string
  annotation: File
  outfilename_stringtie_1: string
  outfilename_stringtie_2: string
  program: File
  script: File
  metadata: File

outputs:
  star_readmap_out:
    type: Directory
    outputSource: star_folder/out
  samtools_out:
    type: Directory
    outputSource: samtools_folder/out
  stringtie_out:
    type: Directory
    outputSource: stringtie_folder/out
  prepDE_out:
    type: Directory
    outputSource: prepDE_folder/out
  DESeq2_out:
    type: Directory
    outputSource: DESeq2_folder/out
steps:
# STAR
  star_readmap_1:
    run: ../../cwl-tools/docker/STAR_readmap.cwl
    in:
      Threads: Threads
      genomeDir: genomeDir
      readFilesIn: readFilesIn_1
      outFileNamePrefix: outFileNamePrefix_1
    out: [sam_output, star_read_out]

  star_readmap_2:
    run: ../../cwl-tools/docker/STAR_readmap.cwl
    in:
      Threads: Threads
      genomeDir: genomeDir
      readFilesIn: readFilesIn_2
      outFileNamePrefix: outFileNamePrefix_2
    out: [sam_output, star_read_out]
  star_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        dir1: Directory
        dir2: Directory
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "star",
            "listing": [inputs.dir1, inputs.dir2]
            } };
          }
    in:
      dir1: star_readmap_1/star_read_out
      dir2: star_readmap_2/star_read_out
    out: [out]
  

# Samtools
  samtools_1:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: star_readmap_1/sam_output
      threads: Threads
      threads2: Threads
      outfilename: outfilename_samtools_1
    out: [samtools_out]

  samtools_2:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: star_readmap_2/sam_output
      threads: Threads
      threads2: Threads
      outfilename: outfilename_samtools_2
    out: [samtools_out]
  samtools_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        file1: File
        file2: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "samtools",
            "listing": [inputs.file1, inputs.file2]
            } };
          }
    in:
      file1: samtools_1/samtools_out
      file2: samtools_2/samtools_out
    out: [out]

#Stringtie
  stringtie_1:
    run: ../../cwl-tools/docker/stringtie.cwl
    in:
      input_bam: samtools_1/samtools_out
      threads: Threads
      annotation: annotation
      outfilename: outfilename_stringtie_1
    out: [stringtie_out]

  stringtie_2:
    run: ../../cwl-tools/docker/stringtie.cwl
    in:
      input_bam: samtools_2/samtools_out
      threads: Threads
      annotation: annotation
      outfilename: outfilename_stringtie_2
    out: [stringtie_out]
  stringtie_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        file1: File
        file2: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "stringtie",
            "listing": [inputs.file1, inputs.file2]
            } };
          }
    in:
      file1: stringtie_1/stringtie_out
      file2: stringtie_2/stringtie_out
    out: [out]
  
  prepDE:
    run: ../../cwl-tools/docker/prepDE.cwl
    in:
     program: program
     gtfs: [stringtie_1/stringtie_out, stringtie_2/stringtie_out]
    out: [gene_output, transcript_output]
  
  prepDE_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        file1: File
        file2: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "prepDE",
            "listing": [inputs.file1, inputs.file2]
            } };
          }
    in:
      file1: prepDE/gene_output
      file2: prepDE/transcript_output
    out: [out]

  DESeq2:
    run: ../../cwl-tools/docker/DESeq2.cwl
    in:
      script: script
      count_matrix: prepDE/gene_output
      metadata: metadata
    out: [DESeq2_out]
  
  DESeq2_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        file1: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "DESeq2",
            "listing": [inputs.file1]
            } };
          }
    in:
      file1: DESeq2/DESeq2_out
    out: [out]
