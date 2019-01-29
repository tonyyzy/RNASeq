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
  outFileNamePrefix_1: string?
  outFileNamePrefix_2: string?
  samfile: File
  outfilename_samtools: string
  annotation: File
  outfilename_stringtie_1: string
  outfilename_stringtie_2: string
  program: File
  input_name: string
  input_name_2: string
  name: string
  script: string
  metadata: File

outputs:
  final_out:
    type: File
    outputSource: prepDE/gene_output

steps:
  star_readmap_1:
    run: STAR_readmap.cwl
    in:
      Threads: Threads
      genomeDir: genomeDir
      readFilesIn: readFilesIn_1
      outFileNamePrefix: outFileNamePrefix_1
    out: [sam_output]

  star_readmap_2:
    run: STAR_readmap.cwl
    in:
      Threads: Threads
      genomeDir: genomeDir
      readFilesIn: readFilesIn_2
      outFileNamePrefix: outFileNamePrefix_2
    out: [sam_output]

  samtools_1:
    run: samtools.cwl
    in:
      samfile: star_readmap_1/sam_output
      threads: Threads
      threads2: Threads
      outfilename: outfilename_samtools
    out: [samtools_out]

  samtools_2:
    run: samtools.cwl
    in:
      samfile: star_readmap_2/sam_output
      threads: Threads
      threads2: Threads
      outfilename: outfilename_samtools
    out: [samtools_out]

  stringtie_1:
    run: stringtie.cwl
    in:
      input_bam: samtools_1/samtools_out
      threads: Threads
      annotation: annotation
      outfilename: outfilename_stringtie_1
    out: [stringtie_out]

  stringtie_2:
    run: stringtie.cwl
    in:
      input_bam: samtools_2/samtools_out
      threads: Threads
      annotation: annotation
      outfilename: outfilename_stringtie_2
    out: [stringtie_out]

  prepDE:
    run: prepDE.cwl
    in:
     program: program
     gtfs: [stringtie_1/stringtie_out, stringtie_2/stringtie_out]
    out: [gene_output]

  DESeq2:
    run: DESeq2.cwl
    in:
      script: script
      count_matrix: prepDE/gene_output
      metadata: metadata
    out: [DESeq2_out]
