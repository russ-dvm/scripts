#!/bin/bash
#RSF, rfrase03@uoguelph.ca
#EQUINE
#estimate depth of coverage and on-target reads.

genome=~/genome/genome.fa
capture=~/equine/2014_11_24/nimblegen/3.padded.merged.1000.capture.targets.bed
sorted=~/equine/2014_11_24/sorted
on_target=~/equine/2014_11_24/qc/on_target

# touch "$on_target"/on.target.reads.txt
# echo -e "group"'\t'"reads" >> "$on_target"/on.target.reads.txt

for file in $sorted/*sorted.bam
do
	full_name="${file#~/equine/2014_11_24/sorted/}"
	group="${full_name%_sorted.bam}"
	echo -n $group...
# 	reads=$(bedtools intersect -bed -abam $file -b $capture | wc -l)
# 	echo -e "$group"'\t'"$reads" >> "$on_target"/on.target.reads.txt
	java -Xmx25g -jar ~/java/picard/picard.jar CollectAlignmentSummaryMetrics METRIC_ACCUMULATION_LEVEL=ALL_READS INPUT=$file OUTPUT=~/equine/2014_11_24/qc/picard_metrics/"$group"_picard_metrics.txt REFERENCE_SEQUENCE=$genome VALIDATION_STRINGENCY=LENIENT
	echo done.
# 	java -jar ~/java/GATK/GenomeAnalysisTK.jar -T DepthOfCoverage -R $genome -I $file -o /media/russ/data/bovine/2015_12_04/depth/"$group" -L $capture -ct 1 -ct 10 -ct 20
done
