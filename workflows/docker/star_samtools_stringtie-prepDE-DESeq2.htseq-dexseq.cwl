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
  prepDE_script: File
  DESeq2_script: File
  htseq_prepare_script: File
  htseq_count_script: File
  dexseq_script: File
  metadata: File

outputs:
  star_readmap_out:
    type: Directory
    outputSource: star_folder/out
  samtools_out:
    type: Directory
    outputSource: samtools_folder/out
  stringtie_out:
    type: Directory
    outputSource: stringtie_folder/out
  prepDE_out:
    type: Directory
    outputSource: prepDE_folder/out
  DESeq2_out:
    type: Directory
    outputSource: DESeq2_folder/out
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
  
  star_readmap_5:
    run: ../../cwl-tools/docker/STAR_readmap.cwl
    in:
      threads: threads
      genomeDir: genomeDir
      readFilesIn: fastq5
      outFileNamePrefix: subject_name5
    out: [sam_output, star_read_out]
  star_readmap_6:
    run: ../../cwl-tools/docker/STAR_readmap.cwl
    in:
      threads: threads
      genomeDir: genomeDir
      readFilesIn: fastq6
      outFileNamePrefix: subject_name6
    out: [sam_output, star_read_out]

  star_folder:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        dir1: Directory
        dir2: Directory
        dir3: Directory
        dir4: Directory
        dir5: Directory
        dir6: Directory
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "star",
            "listing": [inputs.dir1,
                        inputs.dir2,
                        inputs.dir3,
                        inputs.dir4,
                        inputs.dir5,
                        inputs.dir6]
            } };
          }
    in:
      dir1: star_readmap_1/star_read_out
      dir2: star_readmap_2/star_read_out
      dir3: star_readmap_3/star_read_out
      dir4: star_readmap_4/star_read_out
      dir5: star_readmap_5/star_read_out
      dir6: star_readmap_6/star_read_out
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

  samtools_5:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: star_readmap_5/sam_output
      threads: threads
      outfilename:
        source: [subject_name5]
        valueFrom: $(self + ".bam")
    out: [samtools_out]

  samtools_6:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: star_readmap_6/sam_output
      threads: threads
      outfilename:
        source: [subject_name6]
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
        file5: File
        file6: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "samtools",
            "listing": [inputs.file1,
                        inputs.file2,
                        inputs.file3,
                        inputs.file4,
                        inputs.file5,
                        inputs.file6]
            } };
          }
    in:
      file1: samtools_1/samtools_out
      file2: samtools_2/samtools_out
      file3: samtools_3/samtools_out
      file4: samtools_4/samtools_out
      file5: samtools_5/samtools_out
      file6: samtools_6/samtools_out
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

  stringtie_5:
    run: ../../cwl-tools/docker/stringtie.cwl
    in:
      input_bam: samtools_5/samtools_out
      threads: threads
      annotation: annotation
      outfilename:
        source: [subject_name5]
        valueFrom: $(self + ".gtf")
    out: [stringtie_out]

  stringtie_6:
    run: ../../cwl-tools/docker/stringtie.cwl
    in:
      input_bam: samtools_6/samtools_out
      threads: threads
      annotation: annotation
      outfilename:
        source: [subject_name6]
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
        file5: File
        file6: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "stringtie",
            "listing": [inputs.file1,
                        inputs.file2,
                        inputs.file3,
                        inputs.file4,
                        inputs.file5,
                        inputs.file6]
            } };
          }
    in:
      file1: stringtie_1/stringtie_out
      file2: stringtie_2/stringtie_out
      file3: stringtie_3/stringtie_out
      file4: stringtie_4/stringtie_out
      file5: stringtie_5/stringtie_out
      file6: stringtie_6/stringtie_out
    out: [out]
  
  prepDE:
    run: ../../cwl-tools/docker/prepDE.cwl
    in:
     program: prepDE_script
     gtfs: [stringtie_1/stringtie_out,
            stringtie_2/stringtie_out,
            stringtie_3/stringtie_out,
            stringtie_4/stringtie_out,
            stringtie_5/stringtie_out,
            stringtie_6/stringtie_out]
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
        valueFrom: "yes"
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
        valueFrom: "yes"
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

  htseq_count_5:
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
      sam: samtools_5/samtools_out
      outname:
        source: [subject_name5]
        valueFrom: $(self + "_htseq_count.csv")
    out: [output]

  htseq_count_6:
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
      sam: samtools_6/samtools_out
      outname:
        source: [subject_name6]
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
        file5: File
        file6: File
      outputs:
        out: Directory
      expression: |
        ${
          return {"out": {
            "class": "Directory",
            "basename": "htseq_count",
            "listing": [inputs.file1,
                        inputs.file2,
                        inputs.file3,
                        inputs.file4,
                        inputs.file5,
                        inputs.file6]
            } };
          }
    in:
      file1: htseq_count_1/output
      file2: htseq_count_2/output
      file3: htseq_count_3/output
      file4: htseq_count_4/output
      file5: htseq_count_5/output
      file6: htseq_count_6/output
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