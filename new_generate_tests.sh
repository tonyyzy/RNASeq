module purge
module load python/3.7.1 python/2.7.11 star stringtie samtools R miso

########################################################################
############################## STAR INDEX ##############################
########################################################################

# Delete existing STARIndex directory
if [ -d "./tests/STARIndex" ]; then
    rm -r ./tests/STARIndex
fi
if [ -d "./STARIndex" ]; then
    rm -r ./STARIndex
fi

# create directory for star index
mkdir STARIndex

# run star index and store in new directory
star --runMode genomeGenerate --genomeFastaFiles ./tests/test.fa --sjdbGTFfile ./tests/test.gtf --genomeDir ./STARIndex --genomeSAindexNbases 5

# clean
mv ./STARIndex ./tests/
rm ./Log.out

########################################################################
########################### STAR READMAP ###############################
########################################################################

# gz test
gzip -c ./tests/test1.1.fastq > ./tests/test1.1.fastq.gz
gzip -c ./tests/test1.2.fastq > ./tests/test1.2.fastq.gz

############################# TEST 1 ###################################

# Delete existing readmap files
if [ -d "./test1"]; then
    rm -rf ./test1
fi

# create directory for output
mkdir test1

# run read map
star --genomeDir ./tests/STARIndex --readFilesIn ./tests/test1.1.fastq ./tests/test1.2.fastq --outFileNamePrefix ./test1/test1 --outSAMstrandField intronMotif

# clean
mv ./test1/test1Aligned.out.sam ./tests/test1.star.sam
rm -rf test1
tail -n +5 ./tests/test1.star.sam > ./tests/test1.star.tail.sam


############################# TEST 2 ###################################

# Delete existing readmap files
if [ -d "./test2"]; then
    rm -rf ./test2
fi 

# create directory for output
mkdir test2

# run read map
star --genomeDir ./tests/STARIndex --readFilesIn ./tests/test2.1.fastq ./tests/test2.2.fastq --outFileNamePrefix ./test2/test2 --outSAMstrandField intronMotif

# clean
mv ./test2/test2Aligned.out.sam ./tests/test2.star.sam
rm -rf test2
tail -n +5 ./tests/test2.star.sam > ./tests/test2.star.tail.sam

############################# TEST 3 ###################################

# Delete existing readmap files
if [ -d "./test3"]; then
    rm -rf ./test3
fi 

# create directory for output
mkdir test3

# run read map
star --genomeDir ./tests/STARIndex --readFilesIn ./tests/test3.fastq --outFileNamePrefix ./test3/test3 --outSAMstrandField intronMotif

# clean
mv ./test3/test3Aligned.out.sam ./tests/test3.star.sam
rm -rf test3
tail -n +5 ./tests/test3.star.sam > ./tests/test3.star.tail.sam

############################# TEST 4 ###################################

# Delete existing readmap files
if [ -d "./test4"]; then
    rm ./test4
fi 

# create directory for output
mkdir test4

# run read map
star --genomeDir ./tests/STARIndex --readFilesIn ./tests/test4.fastq --outFileNamePrefix ./test4/test4 --outSAMstrandField intronMotif

# clean
mv ./test4/test4Aligned.out.sam ./tests/test4.star.sam
rm -rf test4
tail -n +5 ./tests/test4.star.sam > ./tests/test4.star.tail.sam

########################################################################
##########################  STAR SAMTOOLS ##############################
########################################################################

################################ TEST 1 ################################
samtools sort -o ./tests/test1.star.bam ./tests/test1.star.sam

################################ TEST 2 ################################
samtools sort -o ./tests/test2.star.bam ./tests/test2.star.sam

################################ TEST 3 ################################
samtools sort -o ./tests/test3.star.bam ./tests/test3.star.sam

################################ TEST 4 ################################
samtools sort -o ./tests/test4.star.bam ./tests/test4.star.sam

########################################################################
##############################  STAR STRINGTIE #########################
########################################################################

# Delete existing stringtie directory
if [ -d "./tests/star_stringtie" ]; then
    rm -r ./tests/star_stringtie
fi

# make directory for output
mkdir ./tests/star_stringtie

################################ TEST 1 ################################
stringtie -eB ./tests/test1.star.bam -G ./tests/test.gtf -o ./tests/star_stringtie/test1/test1.gtf
tail -n +3 ./tests/star_stringtie/test1/test1.gtf > ./tests/test1.stringtie.gtf

################################ TEST 2 ################################
stringtie -eB ./tests/test2.star.bam -G ./tests/test.gtf -o ./tests/star_stringtie/test2/test2.gtf

################################ TEST 3 ################################
stringtie -eB ./tests/test3.star.bam -G ./tests/test.gtf -o ./tests/star_stringtie/test3/test3.gtf

################################ TEST 4 ################################
stringtie -eB ./tests/test4.star.bam -G ./tests/test.gtf -o ./tests/star_stringtie/test4/test4.gtf

# clean
rm ./tests/star_stringtie/test*/*.ctab

########################################################################
########################  STAR STRINGTIE PREPDE  #######################
########################################################################

# run PrepDE python script
python2.7 ./scripts/prepDE.py -i tests/star_stringtie/

# clean
mv ./gene_count_matrix.csv ./tests/gene_count_matrix.star_prepde.csv
rm ./transcript_count_matrix.csv

########################################################################
########################  STAR STRINGTIE DESEQ2  #######################
########################################################################

# run deseq2 R script
Rscript ./tests/DESeq2.R --count_matrix ./tests/gene_count_matrix.star_prepde.csv --metadata ./tests/metadata.csv --condition condition

# clean
mv ./untreated-treated_DGE_res.csv ./tests/DGE_res.star_prepde.csv
rm -r ./untreated-treated_norm_count.csv

########################################################################
####################  STAR STRINGTIE DESEQ2 FGSEA  #####################
########################################################################

# test for fgsea
Rscript ./tests/GSEA_Script.R --de_res ./tests/DGE_res.star_prepde.csv --gene_set ./tests/reactome.tsv
mv gsea_res.csv ./tests/gsea_res.star_prepde_deseq2.csv 

########################################################################
##############################  HTSEQ PREPARE ##########################
########################################################################

# test for htseq prepare
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_prepare.py ./tests/test.gtf ./tests/test.htseq.gff

########################################################################
###############################  STAR HTSEQ  ###########################
########################################################################

# Delete existing stringtie directory
if [ -d "./tests/star_dexseq" ]; then
    rm -r ./tests/star_dexseq
fi

# create directory for storage
mkdir ./tests/star_dexseq

################################ TEST 1 ################################
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py -p yes -s no -f bam -r pos ./tests/test.htseq.gff ./tests/test1.star.bam ./tests/star_dexseq/test1_htseq_count.csv

################################ TEST 2 ################################
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py -p yes -s no -f bam -r pos ./tests/test.htseq.gff ./tests/test2.star.bam ./tests/star_dexseq/test2_htseq_count.csv

################################ TEST 3 ################################
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py -p no -s no -f bam -r pos ./tests/test.htseq.gff ./tests/test3.star.bam ./tests/star_dexseq/test3_htseq_count.csv

################################ TEST 4 ################################
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py -p no -s no -f bam -r pos ./tests/test.htseq.gff ./tests/test4.star.bam ./tests/star_dexseq/test4_htseq_count.csv

########################################################################
###########################  STAR HTSEQ DEXSEQ  ########################
########################################################################

# run dexseq R script
Rscript ./tests/DEXSeq.R --count_matrix_dir ./tests/star_dexseq --metadata ./tests/metadata.csv --condition condition

# clean
mv untreated-treated_DEE_results.csv ./tests/DEE_results.star_dexseq.csv
rm untreated-treated_norm_count.csv

########################################################################
##############################  MISO INDEX  ############################
########################################################################
perl scripts/gtf2gff3.pl tests/test.gtf | tee tests/test.star_miso.gff

########################################################################
###########################  STAR MISO MERGE  ##########################
########################################################################

# merge untreated bam
samtools merge -f ./tests/test_untreated.star_miso.bam ./tests/test1.star.bam ./tests/test2.star.bam

# merge treated bam
samtools merge -f ./tests/test_treated.star_miso.bam ./tests/test3.star.bam ./tests/test4.star.bam

########################################################################
############################  STAR MISO RUN  ###########################
########################################################################

# Delete existing stringtie directory
if [ -d "./tests/miso_settings.txt" ]; then
    rm -r ./tests/miso_settings.txt
fi

# Delete existing stringtie directory
if [ -d "./tests/MISOindex" ]; then
    rm -r ./tests/MISOindex
fi

touch tests/miso_settings.txt && echo "[data]" >> tests/miso_settings.txt && echo "filter_results = False" >> tests/miso_settings.txt && echo "min_event_reads = 20" >> tests/miso_settings.txt && echo "strand = fr-unstranded" >> tests/miso_settings.txt && echo "" >> tests/miso_settings.txt && echo "[cluster]" >> tests/miso_settings.txt && echo "cluster_command = qsub" >> tests/miso_settings.txt && echo "" >> tests/miso_settings.txt && echo "[sampler]" >> tests/miso_settings.txt && echo "burn_in = 500" >> tests/miso_settings.txt && echo "lag = 10" >> tests/miso_settings.txt && echo "num_iters = 5000" >> tests/miso_settings.txt && echo "num_processors = 4" >> tests/miso_settings.txt


# gtf indexing
index_gff --index tests/test.star_miso.gff ./tests/MISOindex

############################### UNTREATED ##############################

# Delete existing stringtie directory
if [ -d "./tests/untreated_miso" ]; then
    rm -r ./tests/untreated_miso
fi

# Run samtools to get bam files to the correct format
samtools sort ./tests/test_untreated.star_miso.bam -o ./tests/test_untreated.star_miso.sorted.bam
samtools index ./tests/test_untreated.star_miso.sorted.bam

# run miso
miso --settings-filename tests/miso_settings.txt --run ./tests/MISOindex ./tests/test_untreated.star_miso.sorted.bam --output-dir ./tests/untreated_miso --read-len 101 --paired-end 222.0 0
summarize_miso --summarize ./tests/untreated_miso/ ./tests/untreated_miso/summary

################################ TREATED ###############################

# Delete existing stringtie directory
if [ -d "./tests/treated_miso" ]; then
    rm -r ./tests/treated_miso
fi

# Run samtools to get bam files to the correct format
samtools sort ./tests/test_treated.star_miso.bam -o ./tests/test_treated.star_miso.sorted.bam
samtools index ./tests/test_treated.star_miso.sorted.bam

# run miso
miso --settings-filename tests/miso_settings.txt --run ./tests/MISOindex ./tests/test_treated.star_miso.sorted.bam --output-dir ./tests/treated_miso --read-len 101
summarize_miso --summarize ./tests/treated_miso/ ./tests/treated_miso/summary

########################################################################
##########################  STAR MISO COMPARE  #########################
########################################################################

# run miso compare
compare_miso --compare ./tests/untreated_miso/ ./tests/treated_miso/ ./tests/compare_miso

# clean
mv ./tests/compare_miso/untreated_miso_vs_treated_miso/bayes-factors/untreated_miso_vs_treated_miso.miso_bf ./tests/miso_compare.bf
rm -rf ./tests/compare_miso

########################################################################
############################  STAR FEATURECOUNTS  ######################
########################################################################

# Delete existing stringtie directory
if [ -d "./featurecount" ]; then
    rm -r ./featurecount
fi

# set up directory
mkdir featurecount

# set up
cp ./tests/test1.star.bam ./featurecount/test1.bam
cp ./tests/test2.star.bam ./featurecount/test2.bam
cp ./tests/test3.star.bam ./featurecount/test3.bam
cp ./tests/test4.star.bam ./featurecount/test4.bam

# run featurecount
Rscript ./scripts/featurecount.R --d ./featurecount --g ./tests/test.gtf --metadata ./tests/metadata.csv

# clean
mv gene_count_matrix.csv ./tests/gene_count_matrix.star_featurecounts.csv
rm -rf featurecount

########################################################################
######################  STAR FEATURECOUNTS EDGER  ######################
########################################################################

# run edger
Rscript ./tests/EdgeR.R --condition condition --counts ./tests/gene_count_matrix.star_featurecounts.csv --metadata ./tests/metadata.csv

# clean
mv untreated-treated_DGE_res.csv ./tests/DGE_res.star_featurecounts_edger.csv
rm untreated-treated_norm_count.csv

########################################################################
#############################  SALMON INDEX  ###########################
########################################################################

# Delete existing stringtie directory
if [ -d "./tests/SALMONIndex" ]; then
    rm -r ./tests/SALMONIndex
fi

docker run -v $PWD/tests:/tests -t combinelab/salmon:0.12.0 salmon index -t /tests/test.cdna.fa -i /SALMONIndex --type quasi -p 4
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/SALMONIndex ./tests
docker stop $CONTAINER_ID

########################################################################
#############################  SALMON QUANT  ######################
########################################################################

# Delete existing stringtie directory
if [ -d "./tests/salmon_quant" ]; then
    rm -r ./tests/salmon_quant
fi

# create output directory
mkdir ./tests/salmon_quant

################################ TEST 1 ################################
docker run -v $PWD/tests:/tests -t combinelab/salmon:0.12.0 salmon quant --validateMappings -i /tests/SALMONIndex -o /test1 -p 4 -l A -1 /tests/test1.1.fastq -2 /tests/test1.2.fastq
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test1 ./tests/salmon_quant
docker stop $CONTAINER_ID

################################ TEST 2 ################################
docker run -v $PWD/tests:/tests -t combinelab/salmon:0.12.0 salmon quant --validateMappings -i /tests/SALMONIndex -o /test2 -p 4 -l A -1 /tests/test2.1.fastq -2 /tests/test2.2.fastq
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test2 ./tests/salmon_quant
docker stop $CONTAINER_ID

################################ TEST 3 ################################
docker run -v $PWD/tests:/tests -t combinelab/salmon:0.12.0 salmon quant --validateMappings -i /tests/SALMONIndex -o /test3 -p 4 -l A -r /tests/test3.fastq
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test3 ./tests/salmon_quant
docker stop $CONTAINER_ID

################################ TEST 4 ################################
docker run -v $PWD/tests:/tests -t combinelab/salmon:0.12.0 salmon quant --validateMappings -i /tests/SALMONIndex -o /test4 -p 4 -l A -r /tests/test4.fastq
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test4 ./tests/salmon_quant
docker stop $CONTAINER_ID

########################################################################
#############################  SALMON COUNT  #######################
########################################################################

# run salmon count
Rscript ./scripts/salmon_R_script.R --gtf ./tests/test.gtf --metadata ./tests/metadata.csv --salmon_dir ./tests/salmon_quant

# clean
mv gene_count_matrix.csv ./tests/salmon_gene_count.csv
rm -rf gene_abundance_matrix.csv gene_length_matrix.csv

########################################################################
############################  SALMON DESEQ2  #######################
########################################################################

# run deseq2
Rscript ./tests/DESeq2.R --count_matrix ./tests/salmon_gene_count.csv --metadata ./tests/metadata.csv

# clean 
mv untreated-treated_DGE_res.csv ./tests/DGE_res.salmon_deseq2.csv
rm untreated-treated_norm_count.csv

########################################################################
############################## HISAT2 INDEX  ###########################
########################################################################

# Delete existing STARIndex directory
#if [ -d "./tests/HISAT2Index" ]; then
#    rm -r ./HISAT2Index
#fi
#if [ -d "./HISAT2Index" ]; then
#    rm -r ./STARIndex
#fi

########################################################################
############################## HISAT2 ALIGN  ###########################
########################################################################

################################ TEST 1 ################################
docker run -v $PWD/tests:/tests -t quay.io/biocontainers/hisat2:2.1.0--py27h2d50403_2 hisat2 -q -x /tests/ht2_index/test_index -1 /tests/test1.1.fastq -2 /tests/test1.2.fastq --dta-cufflinks -p 1 -S /test1.hisat2.sam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test1.hisat2.sam ./tests
docker stop $CONTAINER_ID

tail -n +5 ./tests/test1.hisat2.sam > ./tests/test1.hisat2.tail.sam

################################ TEST 2 ################################
docker run -v $PWD/tests:/tests -t quay.io/biocontainers/hisat2:2.1.0--py27h2d50403_2 hisat2 -q -x /tests/ht2_index/test_index -1 /tests/test2.1.fastq -2 /tests/test2.2.fastq --dta-cufflinks -p 1 -S /test2.hisat2.sam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test2.hisat2.sam ./tests
docker stop $CONTAINER_ID

################################ TEST 3 ################################
docker run -v $PWD/tests:/tests -t quay.io/biocontainers/hisat2:2.1.0--py27h2d50403_2 hisat2 -q -x /tests/ht2_index/test_index -U /tests/test3.fastq --dta-cufflinks -p 1 -S /test3.hisat2.sam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test3.hisat2.sam ./tests
docker stop $CONTAINER_ID

################################ TEST 4 ################################
docker run -v $PWD/tests:/tests -t quay.io/biocontainers/hisat2:2.1.0--py27h2d50403_2 hisat2 -q -x /tests/ht2_index/test_index -U /tests/test4.fastq --dta-cufflinks -p 1 -S /test4.hisat2.sam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test4.hisat2.sam ./tests
docker stop $CONTAINER_ID

########################################################################
########################  HISAT2 SAMTOOLS ##############################
########################################################################

################################ TEST 1 ################################
samtools sort -o ./tests/test1.hisat2.bam ./tests/test1.hisat2.sam

################################ TEST 2 ################################
samtools sort -o ./tests/test2.hisat2.bam ./tests/test2.hisat2.sam

################################ TEST 3 ################################
samtools sort -o ./tests/test3.hisat2.bam ./tests/test3.hisat2.sam

################################ TEST 4 ################################
samtools sort -o ./tests/test4.hisat2.bam ./tests/test4.hisat2.sam

########################################################################
########################## HISAT2 CUFFLINKS  ###########################
########################################################################

# Delete existing stringtie directory
if [ -d "./tests/hisat2_cufflinks" ]; then
    rm -r ./tests/hisat2_cufflinks
fi

mkdir ./tests/hisat2_cufflinks

################################ TEST 1 ################################
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest cufflinks -G /tests/test.gtf -p 1 -o /test1_cufflinks ./tests/test1.hisat2.bam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test1_cufflinks ./tests/hisat2_cufflinks
docker stop $CONTAINER_ID

################################ TEST 2 ################################
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest cufflinks -G /tests/test.gtf -p 1 -o /test2_cufflinks ./tests/test2.hisat2.bam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test2_cufflinks ./tests/hisat2_cufflinks
docker stop $CONTAINER_ID

################################ TEST 3 ################################
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest cufflinks -G /tests/test.gtf -p 1 -o /test3_cufflinks ./tests/test3.hisat2.bam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test3_cufflinks ./tests/hisat2_cufflinks
docker stop $CONTAINER_ID

################################ TEST 4 ################################
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest cufflinks -G /tests/test.gtf -p 1 -o /test4_cufflinks ./tests/test4.hisat2.bam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test4_cufflinks ./tests/hisat2_cufflinks
docker stop $CONTAINER_ID

########################################################################
########################## HISAT2 CUFFMERGE  ###########################
########################################################################

# arrange first
rm ./tests/assembly_GTF_list.txt
touch ./tests/assembly_GTF_list.txt && echo "/tests/hisat2_cufflinks/test1_cufflinks/transcripts.gtf" >> ./tests/assembly_GTF_list.txt && echo "/tests/hisat2_cufflinks/test2_cufflinks/transcripts.gtf" >> ./tests/assembly_GTF_list.txt && echo "/tests/hisat2_cufflinks/test3_cufflinks/transcripts.gtf" >> ./tests/assembly_GTF_list.txt && echo "/tests/hisat2_cufflinks/test4_cufflinks/transcripts.gtf" >> ./tests/assembly_GTF_list.txt

# index fa file
samtools faidx ./tests/test.fa

# run cuffmerge
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest cuffmerge -s /tests/test.fa -p 1 -o /merged -g /tests/test.gtf /tests/assembly_GTF_list.txt
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/merged/merged.gtf ./tests
docker stop $CONTAINER_ID

########################################################################
########################## HISAT2 CUFFQUANT  ###########################
########################################################################

# Delete existing stringtie directory
if [ -d "./tests/hisat2_cuffquant" ]; then
    rm -r ./tests/hisat2_cuffquant
fi

# create directory for saving
mkdir ./tests/hisat2_cuffquant

################################ TEST 1 ################################
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest cuffquant -p 1 -o /test1 ./tests/merged.gtf ./tests/test1.hisat2.bam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test1 ./tests/hisat2_cuffquant
docker stop $CONTAINER_ID

################################ TEST 2 ################################
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest cuffquant -p 1 -o /test2 ./tests/merged.gtf ./tests/test2.hisat2.bam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test2 ./tests/hisat2_cuffquant
docker stop $CONTAINER_ID

################################ TEST 3 ################################
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest cuffquant -p 1 -o /test3 ./tests/merged.gtf ./tests/test3.hisat2.bam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test3 ./tests/hisat2_cuffquant
docker stop $CONTAINER_ID

################################ TEST 4 ################################
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest cuffquant -p 1 -o /test4 ./tests/merged.gtf ./tests/test4.hisat2.bam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test4 ./tests/hisat2_cuffquant
docker stop $CONTAINER_ID

########################################################################
########################### HISAT2 CUFFDIFF  ###########################
########################################################################

# run cuffdiff
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest cuffdiff -p 1 -L untreated,treated -o /cuffdiff_out --FDR 1 ./tests/merged.gtf ./tests/hisat2_cuffquant/test1/abundances.cxb,./tests/hisat2_cuffquant/test2/abundances.cxb ./tests/hisat2_cuffquant/test3/abundances.cxb,./tests/hisat2_cuffquant/test4/abundances.cxb
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/cuffdiff_out .
docker stop $CONTAINER_ID

# run python script to fix the output
python ./scripts/cuffdiff_file_sort.py cuffdiff_out/gene_exp.diff

# clean
mv ./cuffdiff_out/cuffdiff_out_DGE_res.csv ./tests/DGE_res.hisat_cuffdiff.csv
rm -rf ./cuffdiff_out

########################################################################
########################## HISAT2 TABLEMAKER  ##########################
########################################################################

# Delete existing stringtie directory
if [ -d "./tests/hisat2_tablemaker" ]; then
    rm -r ./tests/hisat2_tablemaker
fi


# create directory for saving
mkdir ./tests/hisat2_tablemaker

################################ TEST 1 ################################
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest tablemaker -p 1 -q -W -G ./tests/merged.gtf -o /test1 /tests/test1.hisat2.bam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test1 ./tests/hisat2_tablemaker
docker stop $CONTAINER_ID

################################ TEST 2 ################################
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest tablemaker -p 1 -q -W -G ./tests/merged.gtf -o /test2 ./tests/test2.hisat2.bam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test2 ./tests/hisat2_tablemaker
docker stop $CONTAINER_ID

################################ TEST 3 ################################
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest tablemaker -p 1 -q -W -G ./tests/merged.gtf -o /test3 ./tests/test3.hisat2.bam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test3 ./tests/hisat2_tablemaker
docker stop $CONTAINER_ID

################################ TEST 4 ################################
docker run -v $PWD/tests:/tests -t filipejesus/cufflinks:latest tablemaker -p 1 -q -W -G ./tests/merged.gtf -o /test4 ./tests/test4.hisat2.bam
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/test4 ./tests/hisat2_tablemaker
docker stop $CONTAINER_ID

########################################################################
########################### HISAT2 BALLGOWN  ###########################
########################################################################

docker run -v $PWD:/tests -t quay.io/biocontainers/bioconductor-ballgown:2.14.0--r351_0 Rscript /tests/tests/ballgown.R --data_dir /tests/tests/hisat2_tablemaker --metadata /tests/tests/metadata.csv --condition condition
CONTAINER_ID=$(docker ps -alq)
docker cp $CONTAINER_ID:/untreated-treated_DGE_res.csv .
docker stop $CONTAINER_ID

mv untreated-treated_DGE_res.csv ./tests/DGE_res.hisat2_ballgown.csv

########################################################################
#############################  HISAT2 STRINGTIE ########################
########################################################################

# Delete existing stringtie directory
if [ -d "./tests/hisat2_stringtie" ]; then
    rm -r ./tests/hisat2_stringtie
fi

# make directory for output
mkdir ./tests/hisat2_stringtie

################################ TEST 1 ################################
stringtie -eB ./tests/test1.hisat2.bam -G ./tests/test.gtf -o ./tests/hisat2_stringtie/test1/test1.gtf

################################ TEST 2 ################################
stringtie -eB ./tests/test2.hisat2.bam -G ./tests/test.gtf -o ./tests/hisat2_stringtie/test2/test2.gtf

################################ TEST 3 ################################
stringtie -eB ./tests/test3.hisat2.bam -G ./tests/test.gtf -o ./tests/hisat2_stringtie/test3/test3.gtf

################################ TEST 4 ################################
stringtie -eB ./tests/test4.hisat2.bam -G ./tests/test.gtf -o ./tests/hisat2_stringtie/test4/test4.gtf

# clean
rm ./tests/hisat2_stringtie/test*/*.ctab

########################################################################
#######################  HISAT2 STRINGTIE PREPDE  ######################
########################################################################

# run PrepDE python script
python2.7 ./scripts/prepDE.py -i tests/hisat2_stringtie/

# clean
mv ./gene_count_matrix.csv ./tests/gene_count_matrix.hisat2_prepde.csv
rm ./transcript_count_matrix.csv

########################################################################
#######################  HISAT2 STRINGTIE DESEQ2  ######################
########################################################################

# run deseq2 R script
Rscript ./tests/DESeq2.R --count_matrix ./tests/gene_count_matrix.hisat2_prepde.csv --metadata ./tests/metadata.csv --condition condition

# clean
mv ./untreated-treated_DGE_res.csv ./tests/DGE_res.hisat2_prepde.csv
rm -r ./untreated-treated_norm_count.csv

########################################################################
##############################  HISAT2 HTSEQ  ##########################
########################################################################

# Delete existing stringtie directory
if [ -d "./tests/hisat2_dexseq" ]; then
    rm -r ./tests/hisat2_dexseq
fi

# create directory for storage
mkdir ./tests/hisat2_dexseq

################################ TEST 1 ################################
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py -p yes -s no -f bam -r pos ./tests/test.htseq.gff ./tests/test1.hisat2.bam ./tests/hisat2_dexseq/test1_htseq_count.csv

################################ TEST 2 ################################
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py -p yes -s no -f bam -r pos ./tests/test.htseq.gff ./tests/test2.hisat2.bam ./tests/hisat2_dexseq/test2_htseq_count.csv

################################ TEST 3 ################################
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py -p no -s no -f bam -r pos ./tests/test.htseq.gff ./tests/test3.hisat2.bam ./tests/hisat2_dexseq/test3_htseq_count.csv

################################ TEST 4 ################################
python2.7 ./scripts/Basic_DEXSeq_scripts/dexseq_count_modified.py -p no -s no -f bam -r pos ./tests/test.htseq.gff ./tests/test4.hisat2.bam ./tests/hisat2_dexseq/test4_htseq_count.csv

########################################################################
##########################  HISAT2 HTSEQ DEXSEQ  #######################
########################################################################

# run dexseq R script
Rscript ./tests/DEXSeq.R --count_matrix_dir ./tests/hisat2_dexseq --metadata ./tests/metadata.csv --condition condition

# clean
mv untreated-treated_DEE_results.csv ./tests/DEE_results.hisat2_dexseq.csv
rm untreated-treated_norm_count.csv

###### clean all 
