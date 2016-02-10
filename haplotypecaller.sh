#!/bin/bash

#Script to run GATK Haplotype Caller on all samples (aka groups). Each sample will take ~30 minutes.
#R.S.F. 2015/08/18

#do groups of 4 first
for i in ~/equine/2014_11_24/realigned/group3/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group4/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group5/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group6/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group19/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group15/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group16/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group17/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group7/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group8/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group9/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group11/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group12/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group20/*_realigned_recalibrated.bam
do
	base=${i%_realigned_recalibrated.bam}
	java -Xmx28g -jar /usr/local/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/genome/genome.fa -D ~/equine/snp/dbSNP/merged_chr_dbSNP.vcf -ip 100 -L ~/equine/2014_11_24/ref_files/intervals.intervals --genotyping_mode DISCOVERy -pcrModel CONSERVATIVE -ploidy 8 -ERC GVCF -I $i -o $base.g.vcf
	echo .
	echo .
	echo .
	echo Finished that one, ready for the next!
done

#do groups of 5
for filename in ~/equine/2014_11_24/realigned/group1/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group2/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group24/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group14/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group10/*_realigned_recalibrated.bam ~/equine/2014_11_24/realigned/group13/*_realigned_recalibrated.bam
do
	base=${filename%_realigned_recalibrated.bam}
	java -Xmx28g -jar /usr/local/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/genome/genome.fa -D ~/equine/snp/dbSNP/merged_chr_dbSNP.vcf -ip 100 -L ~/equine/2014_11_24/ref_files/intervals.intervals --genotyping_mode DISCOVERy -pcrModel CONSERVATIVE -ploidy 10 -ERC GVCF -I $filename -o $base.g.vcf
	echo .
	echo .
	echo .
	echo Finished that one, ready for the next!
done

#do groups of 3
for filename in ~/equine/2014_11_24/realigned/group18/*_realigned_recalibrated.bam
do
	base=${filename%_realigned_recalibrated.bam}
	java -Xmx28g -jar /usr/local/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/genome/genome.fa -D ~/equine/snp/dbSNP/merged_chr_dbSNP.vcf -ip 100 -L ~/equine/2014_11_24/ref_files/intervals.intervals --genotyping_mode DISCOVERy -pcrModel CONSERVATIVE -ploidy 6 -ERC GVCF -I $filename -o $base.g.vcf
	echo .
	echo .
	echo .
	echo Finished that one, ready for the next!
done

echo .
echo .
echo .
echo .
echo .
echo Job complete.