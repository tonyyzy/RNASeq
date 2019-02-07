#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
requirements:
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  threads: int
  genomeDir: Directory
  annotation: File
  subject_name1: string
  subject_name2: string
  fastq1: File[]
  fastq2: File[]
  htseq_prepare_script: File
  htseq_count_script: File
  dexseq_script: File
  metadata: File

outputs:
  hisat2_align_out:
    type: Directory
    outputSource: hisat2_align_folder/out
  # samtools_out:
  #   type: Directory
  #   outputSource: samtools_folder/out
  htseq_prepare_out:
    type: Directory
    outputSource: htseq_prepare_folder/out
  htseq_count_out:
    type: Directory
    outputSource: htseq_count_folder/out
  dexseq_out:
    type: Directory
    outputSource: dexseq_folder/out

steps:
  hisat2_align_1:
    run: ../../cwl-tools/docker/hisat2_align.cwl
    in:
      threads: threads
      index_directory: genomeDir
      first_pair:
        source: fastq1
        valueFrom: $(self[0])
      second_pair:
        source: fastq1
        valueFrom: $(self[1])
      sam_name:
        source: subject_name1
        valueFrom: $(self + ".sam")
    out: [sam_output, hisat2_align_out]
  
  hisat2_align_2:
    run: ../../cwl-tools/docker/hisat2_align.cwl
    in:
      threads: threads
      index_directory: genomeDir
      first_pair:
        source: fastq2
        valueFrom: $(self[0])
      second_pair:
        source: fastq2
        valueFrom: $(self[1])
      sam_name:
        source: subject_name2
        valueFrom: $(self + ".sam")
    out: [sam_output, hisat2_align_out]
  
  hisat2_align_folder:
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
            "basename": "hisat2",
            "listing": [inputs.dir1, inputs.dir2]
            } };
          }
    in:
      dir1: hisat2_align_1/hisat2_align_out
      dir2: hisat2_align_2/hisat2_align_out
    out: [out]

  htseq_prepare:
    run: ../../cwl-tools/docker/htseq.cwl
    in:
      input_script: htseq_prepare_script
      gtf: annotation
      gff_name:
        source: [annotation]
        valueFrom: $(self.nameroot + ".gff")
    out: [output]
  
  htseq_prepare_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        file: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "htseq_prepare",
            "listing": [inputs.file]
            } };
          }
    in:
      file: htseq_prepare/output
    out: [out]

  htseq_count_1:
    run: ../../cwl-tools/docker/htseq.cwl
    in:
      input_script: htseq_count_script
      gff: htseq_prepare/output
      sam: hisat2_align_1/sam_output
      outname:
        source: [subject_name1]
        valueFrom: $(self + "_htseq_count.csv")
    out: [output]

  htseq_count_2:
    run: ../../cwl-tools/docker/htseq.cwl
    in:
      input_script: htseq_count_script
      gff: htseq_prepare/output
      sam: hisat2_align_2/sam_output
      outname:
        source: [subject_name2]
        valueFrom: $(self + "_htseq_count.csv")
    out: [output]

  htseq_count_folder:
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
            "basename": "htseq_count",
            "listing": [inputs.file1, inputs.file2]
            } };
          }
    in:
      file1: htseq_count_1/output
      file2: htseq_count_2/output
    out: [out]
  
  dexseq:
    run: ../../cwl-tools/docker/dexseq.cwl
    in:
      input_script: dexseq_script
      counts_matrix: htseq_count_folder/out
      gff: htseq_prepare_folder/out
      metadata: metadata
    out: [output]
  
  dexseq_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        file: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "DEXSeq",
            "listing": [inputs.file]
            } };
          }
    in:
      file: dexseq/output
    out: [out]