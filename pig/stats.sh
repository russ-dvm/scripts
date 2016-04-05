#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#estimate depth of coverage and on-target reads.

genome=/media/russ/data/porcine/genome/Sequence/Custom/genome.custom.fa
intervals=/media/russ/data/porcine/ref_files/intervals.intervals
insert=/media/russ/data/porcine/insert_distribution
capture=/media/russ/data/porcine/ref_files/capture_padded_final.bed
sorted=/media/russ/data/porcine/3_aligned
on_target=/media/russ/data/porcine/on_target
picard=/media/russ/data/porcine/picard


touch "$on_target"/"$1"_on.target.reads.txt
echo -e "group"'\t'"reads" >> "$on_target"/"$1"_on.target.reads.txt

for file in "$1"*de-dup.bam
do
	group="${file%_de-dup.bam}"
	echo -n $group...
	reads=$(bedtools intersect -bed -abam $file -b $capture | wc -l)
	echo -e "$group"'\t'"$reads" >> "$on_target"/"$1"_on.target.reads.txt
	echo done.
# 	java -Xmx25g -jar ~/java/GATK/GenomeAnalysisTK.jar -T DepthOfCoverage -R $genome -I $file -o /media/russ/data/porcine/depth/"$group" -L $capture
	java -Xmx25g -jar ~/java/picard/picard.jar CollectAlignmentSummaryMetrics METRIC_ACCUMULATION_LEVEL=SAMPLE INPUT=$file OUTPUT="$picard"/"$group"_picard_metrics.txt REFERENCE_SEQUENCE=$genome VALIDATION_STRINGENCY=LENIENT

done

#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#estimate insert size distribution using picard.

# realigned_dir=/media/russ/data/bovine/2015_12_13/6_realigned
# genome=/media/russ/data/bovine/genome/Sequence/WholeGenomeFasta/genome.chr.fa
# intervals=/media/russ/data/bovine/ref_files/intervals.intervals
# insert=/media/russ/data/bovine/2015_12_13/insert_distribution
# 
# for file in $realigned_dir/*recal.bam
# do
# 	full_name="${file#/media/russ/data/bovine/2015_12_13/6_realigned/}"
# 	group="${full_name%_realigned_recal.bam}"
# 	java -Xmx25g -jar ~/java/picard/picard.jar CollectInsertSizeMetrics HISTOGRAM_FILE="$insert"/"$group"_plot.pdf INPUT=$file OUTPUT="$insert"/"$group"_insert_metrics.txt VALIDATION_STRINGENCY=LENIENT METRIC_ACCUMULATION_LEVEL=SAMPLE
# done