#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#estimate depth of coverage and on-target reads.

genome=/media/russ/data/bovine/genome/Sequence/WholeGenomeFasta/genome.chr.fa
intervals=/media/russ/data/bovine/ref_files/intervals.intervals
dbsnp=/media/russ/data/bovine/dbsnp/intervals.vcf
gvcf=/media/russ/data/bovine/2015_12_13/8_gvcfs
insert=/media/russ/data/bovine/2015_12_13/insert_distribution
capture=/media/russ/data/bovine/ref_files/padded_capture_coord.bed
sorted=/media/russ/data/bovine/2015_12_13/5_sorted
on_target=/media/russ/data/bovine/2015_12_13/on_target

touch "$on_target"/on.target.reads.txt
echo -e "group"'\t'"reads" >> "$on_target"/on.target.reads.txt

for file in $sorted/*de-dup.bam
do
	full_name="${file#/media/russ/data/bovine/2015_12_13/5_sorted/}"
	group="${full_name%_de-dup.bam}"
	echo -n $group...
	reads=$(bedtools intersect -bed -abam $file -b $capture | wc -l)
	echo -e "$group"'\t'"$reads" >> "$on_target"/on.target.reads.txt
	echo done.
	java -jar ~/java/GATK/GenomeAnalysisTK.jar -T DepthOfCoverage -R $genome -I $file -o /media/russ/data/bovine/2015_12_13/depth/"$group" -L $capture -ct 1 -ct 10 -ct 20
done
