# important things to remember
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
      threads: Threads
      threads2: Threads
      outfilename:
        source: [outFileNamePrefix_1]
        valueFrom: |
          ${return(self + ".bam")}
    out: [samtools_out]
```
