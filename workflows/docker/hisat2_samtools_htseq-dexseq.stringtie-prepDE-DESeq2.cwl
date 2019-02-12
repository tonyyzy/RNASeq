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
  prepDE_script: File
  DESeq2_script: File
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
  stringtie_out:
    type: Directory
    outputSource: stringtie_folder/out
  prepDE_out:
    type: Directory
    outputSource: prepDE_folder/out
  DESeq2_out:
    type: Directory
    outputSource: DESeq2_folder/out

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
      sam_name:
        source: subject_name1
        valueFrom: $(self + ".sam")
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
      sam_name:
        source: subject_name2
        valueFrom: $(self + ".sam")
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
      sam_name:
        source: subject_name3
        valueFrom: $(self + ".sam")
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
      sam_name:
        source: subject_name4
        valueFrom: $(self + ".sam")
    out: [sam_output, hisat2_align_out]
  
  hisat2_align_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        dir1: Directory
        dir2: Directory
        dir3: Directory
        dir4: Directory
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "hisat2",
            "listing": [inputs.dir1, inputs.dir2, inputs.dir3, inputs.dir4]
            } };
          }
    in:
      dir1: hisat2_align_1/hisat2_align_out
      dir2: hisat2_align_2/hisat2_align_out
      dir3: hisat2_align_3/hisat2_align_out
      dir4: hisat2_align_4/hisat2_align_out
    out: [out]

# Samtools
  samtools_1:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: hisat2_align_1/sam_output
      threads: threads
      outfilename:
        source: [subject_name1]
        valueFrom: $(self + ".bam")
    out: [samtools_out]

  samtools_2:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: hisat2_align_2/sam_output
      threads: threads
      outfilename:
        source: [subject_name2]
        valueFrom: $(self + ".bam")
    out: [samtools_out]

  samtools_3:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: hisat2_align_3/sam_output
      threads: threads
      outfilename:
        source: [subject_name3]
        valueFrom: $(self + ".bam")
    out: [samtools_out]

  samtools_4:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: hisat2_align_4/sam_output
      threads: threads
      outfilename:
        source: [subject_name4]
        valueFrom: $(self + ".bam")
    out: [samtools_out]

  samtools_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        file1: File
        file2: File
        file3: File
        file4: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "samtools",
            "listing": [inputs.file1, inputs.file2, inputs.file3, inputs.file4]
            } };
          }
    in:
      file1: samtools_1/samtools_out
      file2: samtools_2/samtools_out
      file3: samtools_3/samtools_out
      file4: samtools_4/samtools_out
    out: [out]

  htseq_prepare:
    run: ../../cwl-tools/docker/htseq_prepare.cwl
    in:
      input_script: htseq_prepare_script
      gtf: annotation
      gff_name:
        source: [annotation]
        valueFrom: $(self.nameroot + ".gff")
    out: [output]
  
  htseq_prepare_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        file: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "htseq_prepare",
            "listing": [inputs.file]
            } };
          }
    in:
      file: htseq_prepare/output
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
      sam: samtools_1/samtools_out
      outname:
        source: [subject_name1]
        valueFrom: $(self + "_htseq_count.csv")
    out: [output]


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
      sam: samtools_2/samtools_out
      outname:
        source: [subject_name2]
        valueFrom: $(self + "_htseq_count.csv")
    out: [output]

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
      sam: samtools_3/samtools_out
      outname:
        source: [subject_name3]
        valueFrom: $(self + "_htseq_count.csv")
    out: [output]
  
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
      sam: samtools_4/samtools_out
      outname:
        source: [subject_name4]
        valueFrom: $(self + "_htseq_count.csv")
    out: [output]

  htseq_count_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        file1: File
        file2: File
        file3: File
        file4: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "htseq_count",
            "listing": [inputs.file1, inputs.file2, inputs.file3, inputs.file4]
            } };
          }
    in:
      file1: htseq_count_1/output
      file2: htseq_count_2/output
      file3: htseq_count_3/output
      file4: htseq_count_4/output
    out: [out]
  
  dexseq:
    run: ../../cwl-tools/docker/dexseq.cwl
    in:
      input_script: dexseq_script
      counts_matrix: htseq_count_folder/out
      gff: htseq_prepare_folder/out
      metadata: metadata
    out: [output]
  
  dexseq_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        file: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "DEXSeq",
            "listing": [inputs.file]
            } };
          }
    in:
      file: dexseq/output
    out: [out]

#Stringtie
  stringtie_1:
    run: ../../cwl-tools/docker/stringtie.cwl
    in:
      input_bam: samtools_1/samtools_out
      threads: threads
      annotation: annotation
      outfilename:
        source: [subject_name1]
        valueFrom: $(self + ".gtf")
    out: [stringtie_out]

  stringtie_2:
    run: ../../cwl-tools/docker/stringtie.cwl
    in:
      input_bam: samtools_2/samtools_out
      threads: threads
      annotation: annotation
      outfilename:
        source: [subject_name2]
        valueFrom: $(self + ".gtf")
    out: [stringtie_out]

  stringtie_3:
    run: ../../cwl-tools/docker/stringtie.cwl
    in:
      input_bam: samtools_3/samtools_out
      threads: threads
      annotation: annotation
      outfilename:
        source: [subject_name3]
        valueFrom: $(self + ".gtf")
    out: [stringtie_out]

  stringtie_4:
    run: ../../cwl-tools/docker/stringtie.cwl
    in:
      input_bam: samtools_4/samtools_out
      threads: threads
      annotation: annotation
      outfilename:
        source: [subject_name4]
        valueFrom: $(self + ".gtf")
    out: [stringtie_out]
  stringtie_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        file1: File
        file2: File
        file3: File
        file4: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "stringtie",
            "listing": [inputs.file1, inputs.file2, inputs.file3, inputs.file4]
            } };
          }
    in:
      file1: stringtie_1/stringtie_out
      file2: stringtie_2/stringtie_out
      file3: stringtie_3/stringtie_out
      file4: stringtie_4/stringtie_out
    out: [out]
  
  prepDE:
    run: ../../cwl-tools/docker/prepDE.cwl
    in:
     program: prepDE_script
     gtfs: [stringtie_1/stringtie_out, stringtie_2/stringtie_out, stringtie_3/stringtie_out, stringtie_4/stringtie_out]
    out: [gene_output, transcript_output]
  
  prepDE_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        file1: File
        file2: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "prepDE",
            "listing": [inputs.file1, inputs.file2]
            } };
          }
    in:
      file1: prepDE/gene_output
      file2: prepDE/transcript_output
    out: [out]

  DESeq2:
    run: ../../cwl-tools/docker/DESeq2.cwl
    in:
      script: DESeq2_script
      count_matrix: prepDE/gene_output
      metadata: metadata
    out: [DESeq2_out]
  
  DESeq2_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        file1: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "DESeq2",
            "listing": [inputs.file1]
            } };
          }
    in:
      file1: DESeq2/DESeq2_out
    out: [out]