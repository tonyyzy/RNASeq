module purge
module load python star stringtie samtools

# test for star index
if [ -d "./STARIndex" ]; then
    rm -r ./STARIndex
fi
if [ -d "./tests/STARIndex" ]; then
    rm -r ./tests/STARIndex
fi
mkdir STARIndex
star --runMode genomeGenerate --genomeFastaFiles ./tests/test.fa --sjdbGTFfile ./tests/test.gff3 --sjdbGTFtagExonParentTranscript exon_id --genomeDir ./STARIndex --genomeSAindexNbases 5
mv ./STARIndex ./tests/
rm ./Log.out

# test for star readmap
if [ -d "./test1" ]; then
        rm -r ./test1
fi
mkdir test1
star --genomeDir ./tests/STARIndex --readFilesIn ./tests/test1.1.fastq ./tests/test1.2.fastq --outFileNamePrefix ./test1/test1
cp ./test1/test1Aligned.out.sam ./tests/test1.sam
tail -n +5 ./test1/test1Aligned.out.sam > ./tests/test1.tail.sam
rm -r test1

# test for samtools
samtools view -Su ./tests/test1.sam | samtools sort -o ./tests/test1.bam

# test for stringtie
stringtie ./tests/test1.bam -G ./tests/test.gff3 -o ./test1.stringtie.gtf
tail -n +3 ./test1.stringtie.gtf > ./tests/test1.stringtie.gtf
rm ./test1.stringtie.gtf
