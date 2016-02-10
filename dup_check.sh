#!/bin/bash

#RSF
#2015/11/23
#Ensure the sam files are sorted and then run Picard to mark and remove duplicates. Create a master list of all the metrics output by Picard to then visualize in R. 

touch duplicates_summary.txt
echo -e "Library"'\t'"Unpaired"'\t'"Read_pairs_examined"'\t'"unnmapped_reads"'\t'"unpaired_read_dups"'\t'"read_pair_dups"'\t'"read_pair_optical_dups"'\t'"percent_dup"'\t'"est_library_size" >> duplicates_summary.txt


for file in *.sam
do
	base=${filename%.sam}
	echo Running samtools on $file.	
	samtools sort -O sam -T sam $file -o $base.sorted.sam
	java -jar ~/java/picard/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT=$base.sorted.sam OUTPUT=$base.dup.bam METRICS_FILE=$base_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
	tail -3 $base_metrics.txt | head -1 >> duplicates_summary.txt
done