#!/bin/bash

#RSF, rfrase03@uoguelph.ca, Dec 14

#merge bams from the same group/pool from the two different bovine runs

run1=/media/russ/data/bovine/2015_12_04/8_new_read_groups
run2=/media/russ/data/bovine/2015_12_13/6_realigned
output=/media/russ/data/bovine/merged_runs/1_merged_bams

for bam_run1 in "$run1"/*bam
do
	group="${bam_run1#/media/russ/data/bovine/2015_12_04/8_new_read_groups/}"
	group="${group%.newRG.bam}"

	bam_run2="$run2"/"$group"_realigned_recal.bam


	echo -n Merging $group...
	samtools merge $output/$group.merged.bam $bam_run1 $bam_run2
	echo done.
	echo -n Indexing $group...
	samtools index $output/$group.merged.bam
	echo done.
done