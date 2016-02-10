#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#generate gvcf files on realigned_recalibrated bam files.


realigned_dir=/media/russ/data/bovine/2015_12_04/6_realigned
genome=/media/russ/data/bovine/genome/Sequence/WholeGenomeFasta/genome.chr.fa
intervals=/media/russ/data/bovine/ref_files/intervals.intervals
dbsnp=/media/russ/data/bovine/dbsnp/intervals.vcf
gvcf=/media/russ/data/bovine/2015_12_04/8_gvcfs

for file in $realigned_dir/*.bam
do
	full_name="${file#/media/russ/data/bovine/2015_12_04/6_realigned/}"
	group="${full_name%_realigned_recal.bam}"
	java -Xmx25g -jar ~/java/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $genome -L $intervals -ip 100 -D $dbsnp -I $file --genotyping_mode DISCOVERY -ploidy 5 -stand_call_conf 30 -stand_emit_conf 10 -o "$gvcf"/"$group".g.vcf -ERC GVCF -pcrModel CONSERVATIVE
done
