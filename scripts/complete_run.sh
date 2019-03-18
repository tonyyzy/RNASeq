#!/usr/bin/env sh

# sh complete_run.sh sam exon_size lib_type out_dir

python /usr/local/misopy-0.5.3/misopy/index_gff.py --index $5 indexed_gtf
samtools sort $1 sorted
samtools index sorted.bam

samtools view ./tests/test_untreated.star_miso.sorted.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c | tee log.txt
len=$(awk 'NR == 1 {print $2}' log.txt)
rm log.txt

if [ $3 = "PE" ]
then
  exon_utils --get-const-exons $5 --min-exon-size $2 --output-dir "exon"
  pe_utils --compute-insert-len sorted.bam ./exon/*const_exons.gff --output-dir insert_dist/ | tee log.txt
  grep -A 1 'mean' log.txt > log
  mean=$(awk 'NR == 4 {print $1}' log)
  sd=$(awk 'NR == 4 {print $2}' log)
  rm -rf log.txt log

  miso --run indexed_gtf sorted.bam --output-dir $4 --read-len $len --paired-end $mean $sd
fi
if [ $3 = "SG" ]
then
  miso --run indexed_gtf sorted.bam --output-dir $4 --read-len $len
fi

summarize_miso --summarize $4/ $4/summary

rm -rf $4/batch-logs/*
