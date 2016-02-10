#!/bin/bash

#This script will first create a new file containing only variants that have passed filter for each group.
#Then, it merges these VCF files (pass only) to create normal and diseased groups.
#Better to remove variants that don't pass filter prior to merging, because as long as ONE group has a variant that passes filter, GATK will consider all alleles from all groups.
#By removing failed variants, when the files are merged alleles will be entered as ./././. - better for downstream application in Plink.
#Designed to be run after "depthfilter.sh"
#R.S.F. 2015/07/21

#remove filtered variants from each variant file
echo Getting those passed variants...
for file in group*/*10x.filtered.vcf
do
	base=${file%.10x.filtered.vcf}
	java -jar /usr/local/GATK/GenomeAnalysisTK.jar -T SelectVariants -R ~/genome/genome.fa -V $file -select 'vc.isNotFiltered()' -o $base.10x.pass.only.vcf
done

echo Creating a merged file of all passed variants in normal animals.
java -jar /usr/local/GATK/GenomeAnalysisTK.jar -T CombineVariants -R ~/genome/genome.fa -V:group1 group1/group1.10x.pass.only.vcf -V:group2 group2/group2.10x.pass.only.vcf -V:group3 group3/group3.10x.pass.only.vcf -V:group4 group4/group4.10x.pass.only.vcf -V:group5 group5/group5.10x.pass.only.vcf -V:group6 group6/group6.10x.pass.only.vcf -V:group19 group19/group19.10x.pass.only.vcf -V:group24 group24/group24.10x.pass.only.vcf -o normals.10x.pass.only.vcf

echo Merged file of all passed variants for diseased
java -jar /usr/local/GATK/GenomeAnalysisTK.jar -T CombineVariants -R ~/genome/genome.fa -V:group14 group14/group14.10x.pass.only.vcf -V:group15 group15/group15.10x.pass.only.vcf -V:group16 group16/group16.10x.pass.only.vcf -V:group17 group17/group17.10x.pass.only.vcf -V:group18 group18/group18.10x.pass.only.vcf -V:group7 group7/group7.10x.pass.only.vcf -V:group8 group8/group8.10x.pass.only.vcf -V:group9 group9/group9.10x.pass.only.vcf -V:group10 group10/group10.10x.pass.only.vcf -V:group11 group11/group11.10x.pass.only.vcf -V:group12 group12/group12.10x.pass.only.vcf -V:group13 group13/group13.10x.pass.only.vcf -V:group20 group20/group20.10x.pass.only.vcf -o diseased.10x.pass.only.vcf