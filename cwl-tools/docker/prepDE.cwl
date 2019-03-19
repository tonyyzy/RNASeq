#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: python2

requirements:
  ShellCommandRequirement: {}
  DockerRequirement:
     dockerPull: python:2.7.15-slim

inputs:
  input_script:
    type: File
    inputBinding:
      position: 4
  stringtie_out:
    type: Directory
    inputBinding:
      prefix: -i
      position: 5
  length:
    type: string?
    inputBinding:
      prefix: -l
      position: 6
  cluster:
    type: string?
    inputBinding:
      prefix: -c
      position: 7
  geneID_prefix:
    type: string?
    inputBinding:
      prefix: -s
      position: 8
  cluster_prefix:
    type: string?
    inputBinding:
      prefix: -k
      position: 9
  transcript_2_gene_if_clustering:
    type: string?
    inputBinding:
      prefix: --legend
      position: 10

# if clutering genes that overlap with different gene IDs use "cluster", "cluster_prefix" and "transcript_2_gene_if_clustering" parameters.

outputs:
  gene_count_output:
    type: File
    outputBinding:
      glob: "*gene_count_matrix*"
  transcript_count_output:
    type: File
    outputBinding:
      glob: "*transcript_count_matrix*"
