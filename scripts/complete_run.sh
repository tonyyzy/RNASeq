#!/usr/bin/env sh

# sh complete_run.sh sam exon_size lib_type out_dir

python /usr/local/misopy-0.5.3/misopy/index_gff.py --index $5 indexed_gtf
samtools sort $1 sorted
samtools index sorted.bam
exon_utils --get-const-exons $5 --min-exon-size $2 --output-dir "exon"
pe_utils --compute-insert-len sorted.bam ./exon/*const_exons.gff --output-dir insert_dist/ | tee log.txt

grep -A 1 'mean' log.txt > log

mean=$(awk 'NR == 4 {print $1}' log)
sd=$(awk 'NR == 4 {print $2}' log)
echo $mean
echo $sd

grep -A 1 'min insert' log.txt > log2
len=$(awk 'NR == 1 {print $3}' log2)
echo $len

if [ $3 = "PE" ]
then
  miso --run indexed_gtf sorted.bam --output-dir $4 --read-len $len --paired-end $mean $sd
fi
if [ $3 = "SG" ]
then
  miso --run indexed_gtf sorted.bam --output-dir $4 --read-len $len
fi

summarize_miso --summarize $4/ $4/summary
