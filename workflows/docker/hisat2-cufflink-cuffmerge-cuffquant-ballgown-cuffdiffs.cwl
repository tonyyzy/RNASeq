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
    fasta: File
    subject_name1: string
    subject_name2: string
    subject_name3: string
    subject_name4: string
    fastq1: File[]
    fastq2: File[]
    fastq3: File[]
    fastq4: File[]
    ballgown_script: File
    metadata: File
    conditions: string[]
    Tag: string

outputs:
    hisat2_align_out:
      type: Directory
      outputSource: hisat2_align_folder/out
    samtools_out:
      type: Directory
      outputSource: samtools_folder/out
    cufflink_out:
      type: Directory
      outputSource: cufflink_folder/out
    cuffquant_out:
      type: Directory
      outputSource: cuffquant_folder/out
    cuffmerge_out:
      type: Directory
      outputSource: cuffmerge_folder/out
    cuffdiff_out:
      type: Directory
      outputSource: cuffdiff_folder/out
    tablemaker_out:
      type: Directory
      outputSource: tablemaker_folder/out
    ballgown_out:
      type: Directory
      outputSource: ballgown_folder/out

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
      XSTag: Tag
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
      XSTag: Tag
      output:
        source: subject_name2
        valueFrom: $(self + '.sam')
    out: [sam_output, hisat2_align_out]

  hisat2_align_3:
    run: ../../cwl-tools/docker/hisat2_align.cwl
    in:
      threads: threads
      index_directory: genomeDir
      #first_pair:
      #  source: fastq3
      #  valueFrom: $(self[0])
      #second_pair:
      #  source: fastq3
      #  valueFrom: $(self[1])
      single_file: fastq3
      XSTag: Tag
      output:
        source: subject_name3
        valueFrom: $(self + '.sam')
    out: [sam_output, hisat2_align_out]

  hisat2_align_4:
    run: ../../cwl-tools/docker/hisat2_align.cwl
    in:
      threads: threads
      index_directory: genomeDir
      #first_pair:
      #  source: fastq3
      #  valueFrom: $(self[0])
      #second_pair:
      #  source: fastq3
      #  valueFrom: $(self[1])
      single_file: fastq4
      XSTag: Tag
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

# cufflinks
  cufflinks_1:
    run: ../../cwl-tools/docker/cufflinks.cwl
    in:
      gtf: annotation
      threads: threads
      bam: samtools_1/samtools_out
      output:
        source: [subject_name1]
        valueFrom: $(self + '/')
    out: [cufflink_out, gtf_out]

  cufflinks_2:
    run: ../../cwl-tools/docker/cufflinks.cwl
    in:
      gtf: annotation
      threads: threads
      bam: samtools_2/samtools_out
      output:
        source: [subject_name2]
        valueFrom: $(self + '/')
    out: [cufflink_out, gtf_out]

  cufflinks_3:
    run: ../../cwl-tools/docker/cufflinks.cwl
    in:
      gtf: annotation
      threads: threads
      bam: samtools_3/samtools_out
      output:
        source: [subject_name3]
        valueFrom: $(self + '/')
    out: [cufflink_out, gtf_out]

  cufflinks_4:
    run: ../../cwl-tools/docker/cufflinks.cwl
    in:
      gtf: annotation
      threads: threads
      bam: samtools_4/samtools_out
      output:
        source: [subject_name4]
        valueFrom: $(self + '/')
    out: [cufflink_out, gtf_out]

  cufflink_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item:
        - cufflinks_1/gtf_out
        - cufflinks_2/gtf_out
        - cufflinks_3/gtf_out
        - cufflinks_4/gtf_out
      name:
        valueFrom: "cufflinks"
    out: [out]

# Cuffmerge
  cuffmerge:
    run: ../../cwl-tools/docker/cuffmerge.cwl
    in:
      threads: threads
      gtf: annotation
      fasta: fasta
      cufflinks_output:
        - cufflinks_1/gtf_out
        - cufflinks_2/gtf_out
        - cufflinks_3/gtf_out
        - cufflinks_4/gtf_out
      output:
        valueFrom: "cuffmerge"
    out: [cuffmerge_out, merged_gtf]

  cuffmerge_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item: cuffmerge/merged_gtf
      name:
        valueFrom: "cuffmerge"
    out: [out]

  cuffquant_1:
    run: ../../cwl-tools/docker/cuffquant.cwl
    in:
      threads: threads
      merged_gtf: cuffmerge/merged_gtf
      alignment_file: samtools_1/samtools_out
      output:
        source: [subject_name1]
        valueFrom: $(self)
    out: [cuffquant_out, cxb]

  cuffquant_2:
    run: ../../cwl-tools/docker/cuffquant.cwl
    in:
      threads: threads
      merged_gtf: cuffmerge/merged_gtf
      alignment_file: samtools_2/samtools_out
      output:
        source: [subject_name2]
        valueFrom: $(self)
    out: [cuffquant_out, cxb]

  cuffquant_3:
    run: ../../cwl-tools/docker/cuffquant.cwl
    in:
      threads: threads
      merged_gtf: cuffmerge/merged_gtf
      alignment_file: samtools_3/samtools_out
      output:
        source: [subject_name3]
        valueFrom: $(self)
    out: [cuffquant_out, cxb]

  cuffquant_4:
    run: ../../cwl-tools/docker/cuffquant.cwl
    in:
      threads: threads
      merged_gtf: cuffmerge/merged_gtf
      alignment_file: samtools_4/samtools_out
      output:
        source: [subject_name4]
        valueFrom: $(self)
    out: [cuffquant_out, cxb]

  cuffquant_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item:
        - cuffquant_1/cxb
        - cuffquant_2/cxb
        - cuffquant_3/cxb
        - cuffquant_4/cxb
      name:
        valueFrom: "cuffquant"
    out: [out]

  cuffdiff:
    run: ../../cwl-tools/docker/cuffdiff.cwl
    in:
      threads: threads
      merged_gtf: cuffmerge/merged_gtf
      #libType:
      #  valueFrom: "fr-unstranded"
      #libNorm:
      #  valueFrom: "classic-fpkm"
      FDR:
        valueFrom: "1"
      label: conditions
      condition1_files:
        - cuffquant_1/cxb
        - cuffquant_2/cxb
      condition2_files:
        - cuffquant_3/cxb
        - cuffquant_4/cxb
      output:
        valueFrom: "cuffdiff"
    out: [cuffdiff_out]

  cuffdiff_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item: 
        - cuffdiff/cuffdiff_out
        - cuffdiff/cuffdiff_out
      name:
        valueFrom: "cuffdiff"
    out: [out]

  tablemaker_1:
    run: ../../cwl-tools/docker/tablemaker.cwl
    in:
      threads: threads
      merged_gtf: cuffmerge/merged_gtf
      bam: samtools_1/samtools_out
      output:
        source: [subject_name1]
        valueFrom: $(self)
    out: [tablemaker_out]
  tablemaker_2:
    run: ../../cwl-tools/docker/tablemaker.cwl
    in:
      threads: threads
      merged_gtf: cuffmerge/merged_gtf
      bam: samtools_2/samtools_out
      output:
        source: [subject_name2]
        valueFrom: $(self)
    out: [tablemaker_out]

  tablemaker_3:
    run: ../../cwl-tools/docker/tablemaker.cwl
    in:
      threads: threads
      merged_gtf: cuffmerge/merged_gtf
      bam: samtools_3/samtools_out
      output:
        source: [subject_name3]
        valueFrom: $(self)
    out: [tablemaker_out]

  tablemaker_4:
    run: ../../cwl-tools/docker/tablemaker.cwl
    in:
      threads: threads
      merged_gtf: cuffmerge/merged_gtf
      bam: samtools_4/samtools_out
      output:
        source: [subject_name4]
        valueFrom: $(self)
    out: [tablemaker_out]

  tablemaker_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item:
        - tablemaker_1/tablemaker_out
        - tablemaker_2/tablemaker_out
        - tablemaker_3/tablemaker_out
        - tablemaker_4/tablemaker_out
      name:
        valueFrom: "tablemaker"
    out: [out]

  ballgown:
    run: ../../cwl-tools/docker/ballgown.cwl
    in:
      input_script: ballgown_script
      metadata: metadata
      condition:
        valueFrom: "condition"
      tablemaker_output:
        - tablemaker_folder/out
    out: [gene_matrix, transcript_matrix]

  ballgown_folder:
    run: ../../cwl-tools/folder.cwl
    in:
      item:
        - ballgown/gene_matrix
        - ballgown/transcript_matrix
      name:
        valueFrom: "ballgown"
    out: [out]
