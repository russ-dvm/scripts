#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#calculate basic mapping metrics using picard.

genome=/media/russ/data/bovine/genome/Sequence/WholeGenomeFasta/genome.chr.fa
intervals=/media/russ/data/bovine/ref_files/intervals.intervals
dbsnp=/media/russ/data/bovine/dbsnp/intervals.vcf
picard=/media/russ/data/bovine/2015_12_13/5a_picard
sorted=/media/russ/data/bovine/2015_12_13/5_sorted

for file in $sorted/*de-dup.bam
do
	full_name="${file#/media/russ/data/bovine/2015_12_13/5_sorted/}"
	group="${full_name%_de-dup.bam}"
	java -Xmx25g -jar ~/java/picard/picard.jar CollectAlignmentSummaryMetrics METRIC_ACCUMULATION_LEVEL=ALL_READS INPUT=$file OUTPUT="$picard"/"$group"_picard_metrics.txt REFERENCE_SEQUENCE=$genome VALIDATION_STRINGENCY=LENIENT
done