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
  htseq_prepare_script: File
  htseq_count_script: File
  dexseq_script: File
  metadata: File

outputs:
  hisat2_align_out:
    type: Directory
    outputSource: hisat2_align_folder/out
  samtools_out:
    type: Directory
    outputSource: samtools_folder/out
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
      output:
        source: subject_name1
        valueFrom: $(self + '.sam')
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
      output:
        source: subject_name2
        valueFrom: $(self + '.sam')
    out: [sam_output, hisat2_align_out]

  hisat2_align_3:
    run: ../../cwl-tools/docker/hisat2_align.cwl
    in:
      threads: threads
      index_directory: genomeDir
      # first_pair:
      #   source: fastq2
      #   valueFrom: $(self[0])
      # second_pair:
      #   source: fastq2
      #   valueFrom: $(self[1])
      single_file: fastq3
      output:
        source: subject_name3
        valueFrom: $(self + '.sam')
    out: [sam_output, hisat2_align_out]

  hisat2_align_4:
    run: ../../cwl-tools/docker/hisat2_align.cwl
    in:
      threads: threads
      index_directory: genomeDir
      # first_pair:
      #   source: fastq2
      #   valueFrom: $(self[0])
      # second_pair:
      #   source: fastq2
      #   valueFrom: $(self[1])
      single_file: fastq4
      output:
        source: subject_name4
        valueFrom: $(self + '.sam')
    out: [sam_output, hisat2_align_out]

  hisat2_align_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item:
      - hisat2_align_1/hisat2_align_out
      - hisat2_align_2/hisat2_align_out
      - hisat2_align_3/hisat2_align_out
      - hisat2_align_4/hisat2_align_out
      name:
        valueFrom: "hisat2"
    out: [out]

# Samtools
  samtools_1:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: hisat2_align_1/sam_output
      threads: threads
      outfilename:
        source: [subject_name1]
        valueFrom: $(self + '.bam')
    out: [samtools_out]

  samtools_2:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: hisat2_align_2/sam_output
      threads: threads
      outfilename:
        source: [subject_name2]
        valueFrom: $(self + '.bam')
    out: [samtools_out]

  samtools_3:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: hisat2_align_3/sam_output
      threads: threads
      outfilename:
        source: [subject_name3]
        valueFrom: $(self + '.bam')
    out: [samtools_out]

  samtools_4:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: hisat2_align_4/sam_output
      threads: threads
      outfilename:
        source: [subject_name4]
        valueFrom: $(self + '.bam')
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

  htseq_prepare:
    run: ../../cwl-tools/docker/htseq_prepare.cwl
    in:
      input_script: htseq_prepare_script
      gtf: annotation
      gff_name:
        source: [annotation]
        valueFrom: $(self.nameroot + '.gff')
    out: [output]

  htseq_prepare_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item: htseq_prepare/output
      name:
        valueFrom: "htseq_prepare"
    out: [out]

  htseq_count_1:
    run: ../../cwl-tools/docker/htseq_count.cwl
    in:
      input_script: htseq_count_script
      pairedend:
        valueFrom: "yes"
      stranded:
        valueFrom: "no"
      input_format:
        valueFrom: "bam"
      sorted_by:
        valueFrom: "pos"
      gff: htseq_prepare/output
      bam: samtools_1/samtools_out
      outname:
        source: [subject_name1]
        valueFrom: $(self + '_htseq_count.csv')
    out: [exon_count_output]


  htseq_count_2:
    run: ../../cwl-tools/docker/htseq_count.cwl
    in:
      input_script: htseq_count_script
      pairedend:
        valueFrom: "yes"
      stranded:
        valueFrom: "no"
      input_format:
        valueFrom: "bam"
      sorted_by:
        valueFrom: "pos"
      gff: htseq_prepare/output
      bam: samtools_2/samtools_out
      outname:
        source: [subject_name2]
        valueFrom: $(self + '_htseq_count.csv')
    out: [exon_count_output]

  htseq_count_3:
    run: ../../cwl-tools/docker/htseq_count.cwl
    in:
      input_script: htseq_count_script
      pairedend:
        valueFrom: "no"
      stranded:
        valueFrom: "no"
      input_format:
        valueFrom: "bam"
      sorted_by:
        valueFrom: "pos"
      gff: htseq_prepare/output
      bam: samtools_3/samtools_out
      outname:
        source: [subject_name3]
        valueFrom: $(self + '_htseq_count.csv')
    out: [exon_count_output]

  htseq_count_4:
    run: ../../cwl-tools/docker/htseq_count.cwl
    in:
      input_script: htseq_count_script
      pairedend:
        valueFrom: "no"
      stranded:
        valueFrom: "no"
      input_format:
        valueFrom: "bam"
      sorted_by:
        valueFrom: "pos"
      gff: htseq_prepare/output
      bam: samtools_4/samtools_out
      outname:
        source: [subject_name4]
        valueFrom: $(self + '_htseq_count.csv')
    out: [exon_count_output]

  htseq_count_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item:
      - htseq_count_1/exon_count_output
      - htseq_count_2/exon_count_output
      - htseq_count_3/exon_count_output
      - htseq_count_4/exon_count_output
      name:
        valueFrom: "htseq_count"
    out: [out]

  dexseq:
    run: ../../cwl-tools/docker/dexseq.cwl
    in:
      input_script: dexseq_script
      counts_matrix: htseq_count_folder/out
      gff: htseq_prepare_folder/out
      metadata: metadata
    out: [DEE_out,norm_out]

  dexseq_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item: 
      - dexseq/DEE_out
      - dexseq/norm_out
      name:
        valueFrom: "DEXSeq"
    out: [out]
