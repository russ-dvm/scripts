#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#convert sam files to bams, sort them, and then remove duplicates.


aligned_dir=/media/russ/data/bovine/2015_12_13/4_aligned
sorted_dir=/media/russ/data/bovine/2015_12_13/5_sorted

touch "$sorted_dir"/duplicates_summary.txt
echo -e "Library"'\t'"Unpaired"'\t'"Read_pairs_examined"'\t'"unnmapped_reads"'\t'"unpaired_read_dups"'\t'"read_pair_dups"'\t'"read_pair_optical_dups"'\t'"percent_dup"'\t'"est_library_size" >> "$sorted_dir"/duplicates_summary.txt

for file in $aligned_dir/*
do
	full_name="${file#/media/russ/data/bovine/2015_12_13/4_aligned/}"
	group="${full_name%_unsorted.sam}"
	echo Working on $group...
 	samtools view -huS $file | samtools sort - "$sorted_dir"/"$group"_sorted
 	samtools index "$sorted_dir"/"$group"_sorted.bam
 	java -jar ~/java/picard/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT="$sorted_dir"/"$group"_sorted.bam OUTPUT="$sorted_dir"/"$group"_de-dup.bam METRICS_FILE="$sorted_dir"/metrics/"$group"_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
	head -8 "$sorted_dir"/metrics/"$group"_metrics.txt | tail -1 >> "$sorted_dir"/duplicates_summary.txt
  	samtools index "$sorted_dir"/"$group"_de-dup.bam
 	echo done.
done