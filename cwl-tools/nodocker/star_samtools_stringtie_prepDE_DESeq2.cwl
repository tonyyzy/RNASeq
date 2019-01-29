#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
inputs:
  Threads: string
  genomeDir: Directory
  readFilesIn: File[]
  outFileNamePrefix: string?
  threads: int
  samfile: File
  outfilename_samtools: string
  annotation: File
  outfilename_stringtie: string
  program: string
  input_name: string
  name: string



outputs:
  prepDE_out:
    type: File
    outputSource: prepDE/gene_output

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
      threads: threads
      threads2: threads
      outfilename: outfilename_samtools
    out: [samtools_out]

  stringtie:
    run: stringtie.cwl
    in:
      input_bam: samtools/samtools_out
      threads: threads
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
