# important things to remember
# non-docker version of workflows are deprecated, primary focus is on dockered version
## use inlinejavascript in workflow's step in
1. do `source: [file]` to put files into the step
2. use `self` to refer to the file in `source`
e.g.
```cwl
# Samtools
  samtools_1:
    run: ../../cwl-tools/docker/samtools.cwl
    in:
      samfile: star_readmap_1/sam_output
      threads: threads
      threads2: threads
      outfilename:
        source: [outFileNamePrefix_1]
        valueFrom: |
          ${return(self + ".bam")}
    out: [samtools_out]
```

## Group into directory
```cwl
samtools_folder:
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
          "basename": "samtools",
          "listing": [inputs.file1, inputs.file2]
          } };
        }
  in:
    file1: samtools_1/samtools_out
    file2: samtools_2/samtools_out
  out: [out]
```
