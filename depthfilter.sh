#!/bin/bash

#Simple script to apply a Depth Filter to all *_filtered.vcf files in the directory. 
#The expected directory structure is groupN/groupN_filtered.vcf. The script is designed to be run in the directory that contains all group directories.
#R.S.F. 2015/07/21

#do groups of 4 first
for filename in group3/*.no.depth.vcf group4/*.no.depth.vcf group5/*.no.depth.vcf group6/*.no.depth.vcf group19/*.no.depth.vcf group15/*.no.depth.vcf group16/*.no.depth.vcf group17/*.no.depth.vcf group7/*.no.depth.vcf group8/*.no.depth.vcf group9/*.no.depth.vcf group11/*.no.depth.vcf group12/*.no.depth.vcf group20/*.no.depth.vcf
do
	base=${filename%.no.depth.vcf}
	java -jar /usr/local/GATK/GenomeAnalysisTK.jar -T VariantFiltration -R ~/genome/genome.fa -V $filename --filterExpression "DP < 40" --filterName DepthFilter -o $base.10x.filtered.vcf
	echo GATK complete.
	echo Renaming files...
#	mv ${filename} $base.no.depth.vcf
	echo All finished.
done

#do groups of 5
for filename in group1/*.no.depth.vcf group2/*.no.depth.vcf group24/*.no.depth.vcf group14/*.no.depth.vcf group10/*.no.depth.vcf group13/*.no.depth.vcf
do
	base=${filename%.no.depth.vcf}
	echo Running GATK...
	java -jar /usr/local/GATK/GenomeAnalysisTK.jar -T VariantFiltration -R ~/genome/genome.fa -V $filename --filterExpression "DP < 50" --filterName DepthFilter -o $base.10x.filtered.vcf
	echo GATK complete.
	echo Renaming files...
#	mv ${filename} $base.no.depth.vcf
	echo All finished.
done

#do groups of 3
for filename in group18/*no.depth.vcf
do
	base=${filename%.no.depth.vcf}
	echo Running GATK...
	java -jar /usr/local/GATK/GenomeAnalysisTK.jar -T VariantFiltration -R ~/genome/genome.fa -V $filename --filterExpression "DP < 30" --filterName DepthFilter -o $base.10x.filtered.vcf
	echo GATK complete.
	echo Renaming files...
#	mv ${filename} $base.no.depth.vcf
	echo All finished.
done

echo .
echo .
echo .
echo .
echo .
echo "Job's done."