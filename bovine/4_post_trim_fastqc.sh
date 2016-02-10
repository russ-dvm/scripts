#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#run fastqc on all the fastq files.

trimmed_dir=/media/russ/data/bovine/2015_12_13/3_trimmed
for file in $trimmed_dir/*r1_trimmed.fastq
do
	read1="${file#/media/russ/data/bovine/2015_12_13/3_trimmed/}"
	group_num="${read1%1_trimmed.fastq}"
	read2="$group_num"2_trimmed.fastq
	fastqc --nogroup --noextract "$trimmed_dir"/"$read1" "$trimmed_dir"/"$read2" -o /media/russ/data/bovine/2015_12_13/fastqc/post-trim/
done