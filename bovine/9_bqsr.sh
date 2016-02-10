#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#run gatk's base score recalibration on the realigned bams.


realigned_dir=/media/russ/data/bovine/2015_12_13/6_realigned
bqsr_dir=/media/russ/data/bovine/2015_12_13/7_bqsr
genome=/media/russ/data/bovine/genome/Sequence/WholeGenomeFasta/genome.chr.fa
intervals=/media/russ/data/bovine/ref_files/intervals.intervals
dbsnp=/media/russ/data/bovine/dbsnp/intervals.vcf
list=/media/russ/data/bovine/ref_files/all_bams_run2.list


java -Xmx28g -jar ~/java/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R $genome -L $intervals -knownSites $dbsnp -I $list -o "$bqsr_dir"/all.bams.recal_data.table

java -Xmx28g -jar ~/java/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R $genome -L $intervals -knownSites $dbsnp -I $list -BQSR "$bqsr_dir"/all.bams.recal_data.table -o "$bqsr_dir"/all.bams.post_recal_data.table

java -Xmx28g -jar ~/java/GATK/GenomeAnalysisTK.jar -T AnalyzeCovariants -R $genome -L $intervals -before "$bqsr_dir"/all.bams.recal_data.table -after "$bqsr_dir"/all.bams.post_recal_data.table -plots "$bqsr_dir"/recalibration_plots.pdf
	
	
for file in $realigned_dir/*.bam
do
	full_name="${file#/media/russ/data/bovine/2015_12_13/6_realigned/}"
	group="${full_name%.bam}"
	java -Xmx20g -jar ~/java/GATK/GenomeAnalysisTK.jar -T PrintReads -R $genome -L $intervals -BQSR "$bqsr_dir"/all.bams.recal_data.table -I $file -o "$realigned_dir"/"$group"_recal.bam
done