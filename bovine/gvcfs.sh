#!/bin/bash

#RSF 2015/12/21, rfrase03@uoguelph.ca
#Loop GATK HaplotypeCaller over the 24 bam files. 


realigned_dir=/media/russ/data/bovine/merged_runs/3_realigned
gvcf_dir=/media/russ/data/bovine/merged_runs/4_gvcfs
genome=/media/russ/data/bovine/genome/Sequence/WholeGenomeFasta/genome.chr.fa
intervals=/media/russ/data/bovine/ref_files/intervals.intervals
dbsnp=/media/russ/data/bovine/dbsnp/intervals.vcf



for x in "$realigned_dir"/*bam
do
	group="${x#/media/russ/data/bovine/merged_runs/3_realigned/}"
	group="${group%.merged.realigned.bam}"
	java -Xmx32g -jar ~/java/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R "$genome" -L "$intervals" -ip 100 -D "$dbsnp" --genotyping_mode DISCOVERY -ploidy 10 -stand_call_conf 30 -stand_emit_conf 10 -ERC GVCF -pcrModel CONSERVATIVE -I $x -o "$gvcf_dir"/"$group".g.vcf -log "$gvcf_dir"/"$group".log
done