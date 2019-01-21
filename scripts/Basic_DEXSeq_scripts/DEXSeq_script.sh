#!/usr/bin/bash

# ./DEXSeq_script.sh sam_file_path annotation_file_path  name_for_counts_file metadata_file

echo -ne "\nBaseSpace directory ='$1'\n\n"

if [[ $2 == *"gtf"* ]]
	then
		echo "You entered GTF file, processing to GFF"
		find="gtf"
                replace="gff3"
                new_variable=${2//$find/$replace}
                python dexseq_prepare_annotation.py $2 $new_variable
		echo "done"
	else
		echo "You entered GFF file, no preprocessing needed"
		new_variable=$2
		echo "done"
fi

echo "Generating count table"
for dir in $1; 
do
	echo -ne "\nDir'=$dir'\n\n"
	echo $4
	echo $new_variable
	python dexseq_count.py $new_variable $dir $3
done
echo "Done"

echo "Starting DEXSeq"
./DEXSeq.R --count_matrix_dir $PWD --gff_file_dir $new_variable --metadata $4
