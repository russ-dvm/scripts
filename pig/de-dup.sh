#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#convert sam files to bams, sort them, and then remove duplicates.


aligned_dir=/media/russ/data/porcine/3_aligned
sorted_dir=/media/russ/data/porcine/3_aligned



touch "$sorted_dir"/"$1"_duplicates_summary.txt
echo -e "Library"'\t'"Unpaired"'\t'"Read_pairs_examined"'\t'"unnmapped_reads"'\t'"unpaired_read_dups"'\t'"read_pair_dups"'\t'"read_pair_optical_dups"'\t'"percent_dup"'\t'"est_library_size" >> "$sorted_dir"/"$1"_duplicates_summary.txt

for file in "$1"*
do
	group="${file%_unsorted.sam}"
	echo Working on $group...
 	samtools view -huS $file | samtools sort - > "$sorted_dir"/"$group"_sorted.bam
 	samtools index "$sorted_dir"/"$group"_sorted.bam
 	java -Xmx28g -jar ~/java/picard/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT="$sorted_dir"/"$group"_sorted.bam OUTPUT="$sorted_dir"/"$group"_de-dup.bam METRICS_FILE="$sorted_dir"/metrics/"$group"_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
	head -8 "$sorted_dir"/metrics/"$group"_metrics.txt | tail -1 >> "$sorted_dir"/"$1"_duplicates_summary.txt
  	samtools index "$sorted_dir"/"$group"_de-dup.bam
 	echo done.
done