#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#run fastqc on all the fastq files.

uncompressed_dir=/media/russ/data/bovine/2015_12_13/2_uncompressed

for file in $uncompressed_dir/*r1.fastq
do
	read1="${file#/media/russ/data/bovine/2015_12_13/2_uncompressed/}"
	group_num="${read1%1.fastq}"
	read2="$group_num"2.fastq
	fastqc --nogroup --noextract "$uncompressed_dir"/"$read1" "$uncompressed_dir"/"$read2" -o /media/russ/data/bovine/2015_12_13/fastqc/pre-trim/
done
