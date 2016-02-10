#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#map trimmed fastq files to the bovine genome.


trimmed_dir=/media/russ/data/bovine/2015_12_13/3_trimmed
genome=/media/russ/data/bovine/genome/Sequence/WholeGenomeFasta/genome.chr.fa
aligned_dir=/media/russ/data/bovine/2015_12_13/4_aligned

for file in $trimmed_dir/*r1_trimmed.fastq
do
	read1="${file#/media/russ/data/bovine/2015_12_13/3_trimmed/}"
	group_num="${read1%1_trimmed.fastq}"
	read2="$group_num"2_trimmed.fastq
	num_only="${group_num%r}"
	bwa mem -k 15 -t 8 -M -R "@RG\tID:run2_"$num_only"\tSM:"$num_only"\tPL:ILLUMINA\tPU:000000000-AL3TM\tLB:bov_colecs" $genome "$trimmed_dir"/"$read1" "$trimmed_dir"/"$read2" > "$aligned_dir"/"$num_only"_unsorted.sam
done