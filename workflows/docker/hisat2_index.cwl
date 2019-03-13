#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

inputs:
   threads: int
   fasta: File[]
   ht2base: string
   gtf: File[]
outputs:
   ht_out:
      type: Directory
      outputSource: hisat2_build_folder/out
   log_out:
      type: File
      outputSource: hisat2_build/log
   splice_sites_out:
      type: File
      outputSource: hisat2_ss/splice_sites
   exon_out:
      type: File
      outputSource: hisat2_exon/exons
steps:
   hisat2_ss:
      run: ../../cwl-tools/docker/hisat2_ss.cwl
      in:
         gtf: gtf
      out: [splice_sites]
   hisat2_exon:
      run: ../../cwl-tools/docker/hisat2_exon.cwl
      in:
         gtf: gtf
      out: [exons]
   hisat2_build:
      run: ../../cwl-tools/docker/hisat2_build.cwl
      in:
         threads: threads
         splice_sites: hisat2_ss/splice_sites
         exon: hisat2_exon/exons
         fasta: fasta
         output: ht2base
      out: [ht, log]
   hisat2_build_folder:
      run: ../../cwl-tools/folder.cwl
      in:
         item: hisat2_build/ht
         name:
            valueFrom: "HISAT2Index"
      out: [out]