#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#performs gatk indel realignment on aligned, sorted, and de-duped bams. requires known SNP sets.


sorted_dir=/media/russ/data/bovine/2015_12_13/5_sorted
realigned_dir=/media/russ/data/bovine/2015_12_13/6_realigned
genome=/media/russ/data/bovine/genome/Sequence/WholeGenomeFasta/genome.chr.fa
dbsnp=/media/russ/data/bovine/dbsnp/intervals.vcf
intervals=/media/russ/data/bovine/ref_files/intervals.intervals

for file in $sorted_dir/*de-dup.bam
do
	full_name="${file#/media/russ/data/bovine/2015_12_13/5_sorted/}"
	group="${full_name%_de-dup.bam}"
	java -Xmx25g -jar ~/java/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $genome -L $intervals -ip 100 -I $file -o "$realigned_dir"/"$group".intervals -known $dbsnp
	java -Xmx25g -jar ~/java/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $genome -I $file -targetIntervals "$realigned_dir"/"$group".intervals -o "$realigned_dir"/"$group"_realigned.bam -known $dbsnp
done




