#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#estimate insert size distribution using picard.

realigned_dir=/media/russ/data/bovine/2015_12_13/6_realigned
genome=/media/russ/data/bovine/genome/Sequence/WholeGenomeFasta/genome.chr.fa
intervals=/media/russ/data/bovine/ref_files/intervals.intervals
dbsnp=/media/russ/data/bovine/dbsnp/intervals.vcf
gvcf=/media/russ/data/bovine/2015_12_13/8_gvcfs
insert=/media/russ/data/bovine/2015_12_13/insert_distribution

for file in $realigned_dir/*recal.bam
do
	full_name="${file#/media/russ/data/bovine/2015_12_13/6_realigned/}"
	group="${full_name%_realigned_recal.bam}"
	java -Xmx25g -jar ~/java/picard/picard.jar CollectInsertSizeMetrics HISTOGRAM_FILE="$insert"/"$group"_plot.pdf INPUT=$file OUTPUT="$insert"/"$group"_insert_metrics.txt VALIDATION_STRINGENCY=LENIENT METRIC_ACCUMULATION_LEVEL=SAMPLE
done