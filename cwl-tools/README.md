# Mini documentation for each tools
### Specify the inputs and outputs for each cwl script
- ## salmon_count
   #### inputs:  
   `input_script`: R count script  
   `gtf`: annotation file in gtf format  
   `metadata`: metadata.csv  
   `salmon_dir`: directory with subdirectorys of outputs from `salmon_quant`  
   > maybe change this directory name? like `quant_results`?
   #### outputs:
   `gene_abundance_matrix.csv`  
   `gene_count_matrix.csv` * this one used for DGE analysis e.g. DESeq2  
   `gene_length_matrix.csv`

