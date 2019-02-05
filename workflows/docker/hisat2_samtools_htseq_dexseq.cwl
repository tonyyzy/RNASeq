#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
requirements:
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  Threads: int
  genomeDir: Directory
  readFilesIn_1: File[]
  readFilesIn_2: File[]
  outFileNamePrefix_1: string
  outFileNamePrefix_2: string
  annotation: File
  program: File
  script: File
  metadata: File

outputs:
  hisat_align_out:
    type: Directory
    outputSource: hisat_align_folder/out
  samtools_out:
    type: Directory
    outputSource: samtools_folder/out
  htseq_out:
    type: Directory
    outputSource: htseq_folder/out
  dexseq_out:
    type: Directory
    outputSource: dexseq_folder/out

steps:
  hisat_align_1:
    run: ../../cwl-tools/docker/hisat2_align.cwl
    in:
      Threads: Threads
      genomeDir: genomeDir
      readFilesIn: readFilesIn_1
      outFileNamePrefix: outFileNamePrefix_1
    out: [sam_output, star_read_out]