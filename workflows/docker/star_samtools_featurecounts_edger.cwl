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
  subject_name3: string
  subject_name4: string
  fastq1: File[]
  fastq2: File[]
  fastq3: File[]
  fastq4: File[]
  featurecounts_script: File
  EdgeR_script: File
  metadata: File

outputs:
  star_readmap_out:
    type: Directory
    outputSource: star_folder/out
  samtools_out:
    type: Directory
    outputSource: samtools_folder/out
  featurecounts_out:
    type: Directory
    outputSource: featurecounts_folder/out
  edger_out:
    type: Directory
    outputSource: edger_folder/out

steps:
# STAR
  star_readmap_1:
    run: ../../cwl-tools/docker/STAR_readmap.cwl
    in:
      threads: threads
      genomeDir: genomeDir
      readFilesIn: fastq1
      outFileNamePrefix: subject_name1
    out: [sam_output, star_read_out]

  star_readmap_2:
    run: ../../cwl-tools/docker/STAR_readmap.cwl
    in:
      threads: threads
      genomeDir: genomeDir
      readFilesIn: fastq2
      outFileNamePrefix: subject_name2
    out: [sam_output, star_read_out]
  star_readmap_3:
    run: ../../cwl-tools/docker/STAR_readmap.cwl
    in:
      threads: threads
      genomeDir: genomeDir
      readFilesIn: fastq3
      outFileNamePrefix: subject_name3
    out: [sam_output, star_read_out]
  star_readmap_4:
    run: ../../cwl-tools/docker/STAR_readmap.cwl
    in:
      threads: threads
      genomeDir: genomeDir
      readFilesIn: fastq4
      outFileNamePrefix: subject_name4
    out: [sam_output, star_read_out]

  star_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item:
      - star_readmap_1/star_read_out
      - star_readmap_2/star_read_out
      - star_readmap_3/star_read_out
      - star_readmap_4/star_read_out
      name:
        valueFrom: "star"
    out: [out]


# Samtools
  samtools_1:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: star_readmap_1/sam_output
      threads: threads
      outfilename:
        source: [subject_name1]
        valueFrom: $(self + ".bam")
    out: [samtools_out]

  samtools_2:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: star_readmap_2/sam_output
      threads: threads
      outfilename:
        source: [subject_name2]
        valueFrom: $(self + ".bam")
    out: [samtools_out]

  samtools_3:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: star_readmap_3/sam_output
      threads: threads
      outfilename:
        source: [subject_name3]
        valueFrom: $(self + ".bam")
    out: [samtools_out]

  samtools_4:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: star_readmap_4/sam_output
      threads: threads
      outfilename:
        source: [subject_name4]
        valueFrom: $(self + ".bam")
    out: [samtools_out]

  samtools_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item:
      - samtools_1/samtools_out
      - samtools_2/samtools_out
      - samtools_3/samtools_out
      - samtools_4/samtools_out
      name:
        valueFrom: "samtools"
    out: [out]

  featurecounts:
    run: ../../cwl-tools/docker/featurecounts.cwl
    in:
      input_script: featurecounts_script
      bam_files:
      - samtools_1/samtools_out
      - samtools_2/samtools_out
      - samtools_3/samtools_out
      - samtools_4/samtools_out
      gtf: annotation
      threads: threads
      metadata: metadata
    out: [gene_count_output]

  featurecounts_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item: featurecounts/gene_count_output
      name:
        valueFrom: "featurecounts"
    out: [out]

  edger:
    run: ../../cwl-tools/docker/edger.cwl
    in:
      input_script: EdgeR_script
      count_matrix: featurecounts/gene_count_output
      metadata: metadata
    out: [edger_out]

  edger_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item: edger/edger_out
      name:
        valueFrom: "edger"
    out: [out]
