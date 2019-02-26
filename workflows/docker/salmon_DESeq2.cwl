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
  salmon_index: Directory
  annotation: File
  subject_name1: string
  subject_name2: string
  subject_name3: string
  subject_name4: string
  fastq1: File[]
  fastq2: File[]
  fastq3: File[]
  fastq4: File[]
  salmon_count_script: File
  DESeq2_script: File
  metadata: File

outputs:
  salmon_quant_out:
    type: Directory
    outputSource: salmon_quant_folder/out
  salmon_count_out:
    type: Directory
    outputSource: salmon_count_folder/out
  DESeq2_out:
    type: Directory
    outputSource: DESeq2_folder/out
steps:
  salmon_quant_1:
    run: ../../cwl-tools/docker/salmon_quant.cwl
    in:
      index_directory: salmon_index
      output: subject_name1
      threads: threads
      first_end_fastq:
        source: [fastq1]
        valueFrom: $(self[0])
      second_end_fastq:
        source: [fastq1]
        valueFrom: $(self[1])
    out: [output]

  salmon_quant_2:
    run: ../../cwl-tools/docker/salmon_quant.cwl
    in:
      index_directory: salmon_index
      output: subject_name2
      threads: threads
      first_end_fastq:
        source: [fastq2]
        valueFrom: $(self[0])
      second_end_fastq:
        source: [fastq2]
        valueFrom: $(self[1])
    out: [output]

  salmon_quant_3:
    run: ../../cwl-tools/docker/salmon_quant.cwl
    in:
      index_directory: salmon_index
      output: subject_name3
      threads: threads
      single_fastq: fastq3
    out: [output]

  salmon_quant_4:
    run: ../../cwl-tools/docker/salmon_quant.cwl
    in:
      index_directory: salmon_index
      output: subject_name4
      threads: threads
      single_fastq: fastq4
    out: [output]

  salmon_quant_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item:
      - salmon_quant_1/output
      - salmon_quant_2/output
      - salmon_quant_3/output
      - salmon_quant_4/output
      name:
        valueFrom: "salmon_quant"
    out: [out]
  
  salmon_count:
    run: ../../cwl-tools/docker/salmon_count.cwl
    in:
      input_script: salmon_count_script
      gtf: annotation
      metadata: metadata
      quant_results: salmon_quant_folder/out
    out: [count, length, abundance]

  salmon_count_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item:
      - salmon_count/count
      - salmon_count/length
      - salmon_count/abundance
      name:
        valueFrom: "salmon_count"
    out: [out]

  DESeq2:
    run: ../../cwl-tools/docker/DESeq2.cwl
    in:
      script: DESeq2_script
      count_matrix: salmon_count/count
      metadata: metadata
    out: [DESeq2_out]
  
  DESeq2_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item: DESeq2/DESeq2_out
      name:
        valueFrom: "DESeq2"
    out: [out]
