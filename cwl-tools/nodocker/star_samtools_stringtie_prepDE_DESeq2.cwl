#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

inputs:
  Threads: int
  genomeDir: Directory
  readFilesIn: File[]
  outFileNamePrefix: string?
  samfile: File
  outfilename_samtools: string
  annotation: File
  outfilename_stringtie: string
  program: string
  input_name: string
  name: string
  script: string
  metadata: File

outputs:
  final_out:
    type: File
    outputSource: DESeq2/DESeq2_out

steps:
  star_readmap:
    run: STAR_readmap.cwl
    in:
      Threads: Threads
      genomeDir: genomeDir
      readFilesIn: readFilesIn
      outFileNamePrefix: outFileNamePrefix
    out: [sam_output]

  samtools:
    run: samtools.cwl
    in:
      samfile: star_readmap/sam_output
      threads: Threads
      threads2: Threads
      outfilename: outfilename_samtools
    out: [samtools_out]

  stringtie:
    run: stringtie.cwl
    in:
      input_bam: samtools/samtools_out
      threads: Threads
      annotation: annotation
      outfilename: outfilename_stringtie
    out: [stringtie_out]

  foldering:
    run: sort_files.cwl
    in:
      files: stringtie/stringtie_out
      input_name: input_name
    out: [foldering_out]

  parenting:
    run: parent_generator.cwl
    in:
      sub_directory: foldering/foldering_out
      name: name
    out: [parenting_out]

  prepDE:
    run: prepDE.cwl
    in:
     program: program
     gtfDir: parenting/parenting_out
    out: [gene_output]

  DESeq2:
    run: DESeq2.cwl
    in:
      script: script
      count_matrix: prepDE/gene_output
      metadata: metadata
    out: [DESeq2_out]
