module purge
module load python/3.7.1 python/2.7.11 star stringtie samtools R

# test for star index
if [ -d "./STARIndex" ]; then
    rm -r ./STARIndex
fi
if [ -d "./tests/GenomeIndex" ]; then
    rm -r ./tests/GenomeIndex
fi
mkdir STARIndex
# star --runMode genomeGenerate --genomeFastaFiles ./tests/test.fa --sjdbGTFfile ./tests/test.gtf --sjdbGTFtagExonParentTranscript exon_id --genomeDir ./STARIndex --genomeSAindexNbases 5
star --runMode genomeGenerate --genomeFastaFiles ./tests/test.fa --sjdbGTFfile ./tests/test.gtf --genomeDir ./STARIndex --genomeSAindexNbases 5
mv ./STARIndex ./tests/GenomeIndex
rm ./Log.out

# test for star readmap
# gz test
gzip -c ./tests/test1.1.fastq > ./tests/test1.1.fastq.gz
gzip -c ./tests/test1.2.fastq > ./tests/test1.2.fastq.gz
# test1
if [ -d "./test1" ]; then
        rm -r ./test1
fi
mkdir test1
star --genomeDir ./tests/GenomeIndex --readFilesIn ./tests/test1.1.fastq ./tests/test1.2.fastq --outFileNamePrefix ./test1/test1 --outSAMstrandField intronMotif
cp ./test1/test1Aligned.out.sam ./tests/test1.sam
tail -n +5 ./test1/test1Aligned.out.sam > ./tests/test1.tail.sam
# test2
if [ -d "./test2" ]; then
            rm -r ./test2
fi
mkdir test2
star --genomeDir ./tests/GenomeIndex --readFilesIn ./tests/test2.1.fastq ./tests/test2.2.fastq --outFileNamePrefix ./test2/test2 --outSAMstrandField intronMotif
cp ./test2/test2Aligned.out.sam ./tests/test2.sam
# test3
if [ -d "./test3" ]; then
    rm -r ./test3
fi
mkdir test3
star --genomeDir ./tests/GenomeIndex --readFilesIn ./tests/test3.fastq --outFileNamePrefix ./test3/test3 --outSAMstrandField intronMotif
cp ./test3/test3Aligned.out.sam ./tests/test3.sam
tail -n +5 ./test3/test3Aligned.out.sam > ./tests/test3.tail.sam
# test4
if [ -d "./test4" ]; then
        rm -r ./test4
fi
mkdir test4
star --genomeDir ./tests/GenomeIndex --readFilesIn ./tests/test4.fastq --outFileNamePrefix ./test4/test4 --outSAMstrandField intronMotif
cp ./test4/test4Aligned.out.sam ./tests/test4.sam
# cleanup
rm -r test1 test2 test3 test4

# test for samtools
samtools sort -o ./tests/test1.bam ./tests/test1.sam
samtools sort -o ./tests/test2.bam ./tests/test2.sam
samtools sort -o ./tests/test3.bam ./tests/test3.sam
samtools sort -o ./tests/test4.bam ./tests/test4.sam

# test for stringtie
if [ -d "./tests/stringtie" ]; then
    rm -r ./tests/stringtie
fi
mkdir ./tests/stringtie
stringtie -eB ./tests/test1.bam -G ./tests/test.gtf -o ./tests/stringtie/test1/test1.gtf
tail -n +3 ./tests/stringtie/test1/test1.gtf > ./tests/test1.stringtie.gtf
stringtie -eB ./tests/test2.bam -G ./tests/test.gtf -o ./tests/stringtie/test2/test2.gtf
stringtie -eB ./tests/test3.bam -G ./tests/test.gtf -o ./tests/stringtie/test3/test3.gtf
stringtie -eB ./tests/test4.bam -G ./tests/test.gtf -o ./tests/stringtie/test4/test4.gtf
rm ./tests/stringtie/test*/*.ctab

# test for prepDE
# cwl-runner --outdir=./tests ./cwl-tools/nodocker/prepDE.cwl ./tests/prepDE.yml
python2.7 ./scripts/prepDE.py -i tests/stringtie/
mv ./gene_count_matrix.csv ./tests/
mv ./transcript_count_matrix.csv ./tests/

# test for DESeq2
cwl-runner --outdir=./test_DESeq2 ./cwl-tools/docker/DESeq2.cwl ./tests/DESeq2.yml
cp ./test_DESeq2/DGE_results.csv ./tests/
rm -r ./test_DESeq2

# test for htseq prepare
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_prepare.py ./tests/test.gtf ./tests/test.gff

# test for htseq counts
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py -p yes -s no -f bam -r pos ./tests/test.gff ./tests/test1.bam ./tests/test1_htseq_count.csv
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py -p yes -s no -f bam -r pos ./tests/test.gff ./tests/test2.bam ./tests/test2_htseq_count.csv
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py -p no -s no -f bam -r pos ./tests/test.gff ./tests/test3.bam ./tests/test3_htseq_count.csv
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py -p no -s no -f bam -r pos ./tests/test.gff ./tests/test4.bam ./tests/test4_htseq_count.csv
# test for dexseq
Rscript ./tests/DEXSeq.R --count_matrix_dir ./tests --gff_file_dir ./tests --metadata ./tests/metadata.csv
mv DEE_results.csv ./tests/test_DEE_results.csv

# test for fgsea
Rscript ./tests/GSEA_Script.R --de_res ./tests/test_dge_results.csv --gene_set ./tests/reactome.tsv --doc_name ./tests/test_gsea_res.csv

# test for hisat_align
cwl-runner --outdir=./test_hisat_align cwl-tools/docker/hisat2_align.cwl tests/hisat2_align.yml
tail -n +5 ./test_hisat_align/test1/test1.sam > ./tests/test1.hisat2.tail.sam
rm -r ./test_hisat_align

# Salmon index
cwl-runner ./cwl-tools/docker/salmon_index.cwl ./tests/salmon_index.yml
mv Salmonindex ./tests/

# Salmon quant
mkdir ./tests/salmon_quant
cwl-runner ./cwl-tools/docker/salmon_quant.cwl ./tests/salmon_quant.yml
tail -n +2 ./test2/quant.sf | awk 'BEGIN{OFS=FS="\t"}{$3=sprintf("%3.0f",$3);$4=sprintf("%3.0f",$4)}1' > ./tests/salmon_quant/test2.sf
cwl-runner ./cwl-tools/docker/salmon_quant.cwl ./tests/salmon_quant.single.yml
tail -n +2 ./test3/quant.sf | awk 'BEGIN{OFS=FS="\t"}{$3=sprintf("%3.0f",$3);$4=sprintf("%3.0f",$4)}1' > ./tests/salmon_quant/test3.sf
rm -r test2 test3
# workflow 4
cwl-runner --outdir=./workflow4 ./workflows/docker/salmon_DESeq2.cwl ./tests/salmon_DESeq2.yml
cp ./workflow4/DESeq2/DGE_results.csv ./tests/salmon_DGE_results.csv

# Salmon count
cp ./workflow4/salmon_count/gene_count_matrix.csv ./tests/salmon_gene_count.csv
rm -rf ./workflow4

# workflow 2
cwl-runner --outdir=./workflow2 ./workflows/docker/hisat2_htseq_dexseq.cwl ./tests/hisat2_htseq_dexseq.yml
cp ./workflow2/DEXSeq/DEE_results.csv ./tests/DEE_results.csv
rm -rf ./workflow2
