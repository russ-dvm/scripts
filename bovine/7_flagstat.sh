#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#run samtools flagstat on all aligned bams, that have gone through picard's pcr duplicate removal function.


sorted_dir=/media/russ/data/bovine/2015_12_13/5_sorted


for file in $sorted_dir/*de-dup.bam
do
	full_name="${file#/media/russ/data/bovine/2015_12_13/5_sorted/}"
	group="${full_name%_de-dup.bam}"
	echo -n Generating metrics on "$group"...
	samtools flagstat $file > /media/russ/data/bovine/2015_12_13/flagstat/"$group"_flagstat.txt
	echo done.
done
echo Script complete.