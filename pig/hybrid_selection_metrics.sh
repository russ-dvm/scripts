#!/bin/bash


genome=/media/russ/data/porcine/genome/Sequence/Custom/genome.custom.fa
bait_bed=/media/russ/data/porcine/ref_files/capture_padded_final.bed
target_bed=/media/russ/data/porcine/ref_files/primary_padded_final.bed
hybrid=/media/russ/data/porcine/hybrid_selection

for file in *bam
do
	base="${file%_rea*}"
	echo $file
	echo $base
	samtools view -H "$file" > "$hybrid"/"$base"_temp.txt
	cat "$target_bed" | gawk '{print $1 "\t" $2+1 "\t" $3 "\t+\tinterval_" NR}' > "$hybrid"/"$base"_target_temp_body.txt
	cat "$hybrid"/"$base"_temp.txt "$hybrid"/"$base"_target_temp_body.txt > "$hybrid"/"$base"_target_intervals.txt
	target="$hybrid"/"$base"_target_intervals.txt
	
	cat "$bait_bed" | gawk '{print $1 "\t" $2+1 "\t" $3 "\t+\tinterval_" NR}' > "$hybrid"/"$base"_bait_temp_body.txt
	cat "$hybrid"/"$base"_temp.txt "$hybrid"/"$base"_bait_temp_body.txt > "$hybrid"/"$base"_bait_intervals.txt
	bait="$hybrid"/"$base"_bait_intervals.txt

	
	java -Xmx28g -jar ~/java/picard/picard.jar CalculateHsMetrics BAIT_INTERVALS="$bait" TARGET_INTERVALS="$target" INPUT="$file" OUTPUT="$hybrid"/"$base"_hs_metrics METRIC_ACCUMULATION_LEVEL=LIBRARY REFERENCE_SEQUENCE="$genome" VALIDATION_STRINGENCY=LENIENT
	
done