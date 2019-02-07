# RNASeq

[![Build Status](https://travis-ci.com/tonyyzy/RNASeq.svg?token=5Fwptxoz1iaezXoMzRSd&branch=master)](https://travis-ci.com/tonyyzy/RNASeq)

# We will probably drop nodocker support soon...it's becoming harder to account for both versions
Death to no docker !!!!!!!!!!

## Notes
DESeq2 requires metadata sample name to match gene_count_matrix sample name
will handle this in workflow to make things consistent

## Some libraries needed to setup local dependencies
```bash
sudo apt-get install zlib1g-dev
sudo apt-get install libncurses5-dev
sudo apt-get install libbz2-dev
sudo apt-get install liblzma-dev
```

## Testing

We use `python` and `pytest` as the testing framework. All the tests are stored in `./tests`. Before committing a new test, you should test your test locally first.

Assuming you have all dependencies met

run `$ pytest ./tests`

You should see something similar to the following:

```
$ pytest tests
============================= test session starts ==============================
platform linux -- Python 3.6.3, pytest-3.3.0, py-1.5.2, pluggy-0.6.0
rootdir: /home/travis/build/tonyyzy/RNASeq, inifile:
collected 1 item
tests/test_STAR_index_nodocker.py .                                      [100%]
=========================== 1 passed in 2.01 seconds ===========================
The command "pytest tests" exited with 0.
```
## Webportal structure  
<pre>
RNASeq  
|  
|-- Genome_Index  
|               |-- H. Sapiens  
|               |             |-- fasta  
|               |             |-- annotation  
|               |             |-- STAR_index  
|               |             |-- HISAT2_index  
|               |              
|               |-- P. Yoelii  
|-- data  
|       |  
|       |-- Session ID  
|                     |-- workflow1  
|                     |            |-- cwl file  
|                     |            |-- yml file
|                     |            |-- analysis result directories
|                     |-- fastq  
|                     |        |-- file1.fastq  
|                     |        |-- file2.fastq  
|                     |  
|                     |-- genome  
|                              |-- fasta  
|                              |-- annotation  
|  
|-- RNASeq (repository)  
</pre>
