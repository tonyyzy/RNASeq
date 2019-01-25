module purge
module load python star stringtie samtools

# test for star index
if [ -d "./STARIndex" ]; then
    rm -r ./STARIndex
fi
if [ -d "./tests/GenomeIndex" ]; then
    rm -r ./tests/GenomeIndex
fi
mkdir STARIndex
star --runMode genomeGenerate --genomeFastaFiles ./tests/test.fa --sjdbGTFfile ./tests/test.gff3 --sjdbGTFtagExonParentTranscript exon_id --genomeDir ./STARIndex --genomeSAindexNbases 5
mv ./STARIndex ./tests/GenomeIndex
rm ./Log.out

# test for star readmap
# test1
if [ -d "./test1" ]; then
        rm -r ./test1
fi
mkdir test1
star --genomeDir ./tests/GenomeIndex --readFilesIn ./tests/test1.1.fastq ./tests/test1.2.fastq --outFileNamePrefix ./test1/test1
cp ./test1/test1Aligned.out.sam ./tests/test1.sam
tail -n +5 ./test1/test1Aligned.out.sam > ./tests/test1.tail.sam
# test2
if [ -d "./test2" ]; then
            rm -r ./test2
fi
mkdir test2
star --genomeDir ./tests/GenomeIndex --readFilesIn ./tests/test2.1.fastq ./tests/test2.2.fastq --outFileNamePrefix ./test2/test2
cp ./test2/test2Aligned.out.sam ./tests/test2.sam
# test3
if [ -d "./test3" ]; then
    rm -r ./test3
fi
mkdir test3
star --genomeDir ./tests/GenomeIndex --readFilesIn ./tests/test3.fastq --outFileNamePrefix ./test3/test3
cp ./test3/test3Aligned.out.sam ./tests/test3.sam
# test4
if [ -d "./test4" ]; then
        rm -r ./test4
fi
mkdir test4
star --genomeDir ./tests/GenomeIndex --readFilesIn ./tests/test4.fastq --outFileNamePrefix ./test4/test4
cp ./test4/test4Aligned.out.sam ./tests/test4.sam
# cleanup
rm -r test1 test2 test3 test4

# test for samtools
samtools view -Su ./tests/test1.sam | samtools sort -o ./tests/test1.bam

# test for stringtie
stringtie ./tests/test1.bam -G ./tests/test.gff3 -o ./test1.stringtie.gtf
tail -n +3 ./test1.stringtie.gtf > ./tests/test1.stringtie.gtf
rm ./test1.stringtie.gtf
