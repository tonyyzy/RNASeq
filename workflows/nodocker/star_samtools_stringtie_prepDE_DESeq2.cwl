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
  star_readmap_1_out:
    type: Directory
    outputSource: star_readmap_1/star_read_out
  star_readmap_2_out:
    type: Directory
    outputSource: star_readmap_2/star_read_out
  samtools_1_out:
    type: File
    outputSource: samtools_1/samtools_out
  samtools_2_out:
    type: File
    outputSource: samtools_2/samtools_out
  stringtie_1_out:
    type: File
    outputSource: stringtie_1/stringtie_out
  stringtie_2_out:
    type: File
    outputSource: stringtie_2/stringtie_out
  prepDE_out:
    type: File[]
    outputSource:
      - prepDE/gene_output
      - prepDE/transcript_output
  DESeq2_out:
    type: File
    outputSource: DESeq2/DESeq2_out


steps:
  star_readmap_1:
    run: ../../cwl-tools/nodocker/STAR_readmap.cwl
    in:
      Threads: Threads
      genomeDir: genomeDir
      readFilesIn: readFilesIn_1
      outFileNamePrefix: outFileNamePrefix_1
    out: [sam_output, star_read_out]

  star_readmap_2:
    run: ../../cwl-tools/nodocker/STAR_readmap.cwl
    in:
      Threads: Threads
      genomeDir: genomeDir
      readFilesIn: readFilesIn_2
      outFileNamePrefix: outFileNamePrefix_2
    out: [sam_output, star_read_out]

  samtools_1:
    run: ../../cwl-tools/nodocker/samtools.cwl
    in:
      samfile: star_readmap_1/sam_output
      threads: Threads
      threads2: Threads
      outfilename: outfilename_samtools_1
    out: [samtools_out]

  samtools_2:
    run: ../../cwl-tools/nodocker/samtools.cwl
    in:
      samfile: star_readmap_2/sam_output
      threads: Threads
      threads2: Threads
      outfilename: outfilename_samtools_2
    out: [samtools_out]

  stringtie_1:
    run: ../../cwl-tools/nodocker/stringtie.cwl
    in:
      input_bam: samtools_1/samtools_out
      threads: Threads
      annotation: annotation
      outfilename: outfilename_stringtie_1
    out: [stringtie_out]

  stringtie_2:
    run: ../../cwl-tools/nodocker/stringtie.cwl
    in:
      input_bam: samtools_2/samtools_out
      threads: Threads
      annotation: annotation
      outfilename: outfilename_stringtie_2
    out: [stringtie_out]

  prepDE:
    run: ../../cwl-tools/nodocker/prepDE.cwl
    in:
     program: program
     gtfs: [stringtie_1/stringtie_out, stringtie_2/stringtie_out]
    out: [gene_output, transcript_output]

  DESeq2:
    run: ../../cwl-tools/docker/DESeq2.cwl
    in:
      script: script
      count_matrix: prepDE/gene_output
      metadata: metadata
    out: [DESeq2_out]
