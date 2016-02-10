#!/bin/bash

#RSF 2015/10/08
#Set of commands to select variants from collagenous lectins only; to separate them out by region (upstream, downstream, coding, introns).

#annotate with snpeff
echo -n Running snpEff...
java -Xmx10G -jar ~/java/snpEff/snpEff.jar eff -i vcf -o vcf EquCab2.78 -interval ~/equine/2014_11_24/ref_files/intervals.txt -stats $1.html $1.vcf > $1.snpeff.vcf
java -Xmx10G -jar ~/java/snpEff/snpEff.jar eff -i vcf -o vcf EquCab2.78 -ud 50000 -interval ~/equine/2014_11_24/ref_files/intervals.txt -stats $1.html $1.vcf > $1.snpeff.50kb.vcf
echo done.

#Select variants from collagenous lectins only -ef = exclude filtered, -env = exclude non variant
echo -n Generating VCFs with variants that are in lectins and have passed...

java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T SelectVariants -R ~/genome/genome.fa -L ~/equine/2014_11_24/ref_files/lectin_intervals.intervals -ef -env -V $1.snpeff.vcf -o $1.lectins.vcf
java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T SelectVariants -R ~/genome/genome.fa -L ~/equine/2014_11_24/ref_files/lectin_intervals.intervals -ef -env -V $1.snpeff.50kb.vcf -o $1.lectins.50kb.vcf

echo done.

echo .
echo .
echo .

#Split them into normals and diseased groups.
echo Splitting into normals and diseased groups...
echo .
echo .
echo .
#Normals
java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T SelectVariants -R ~/genome/genome.fa -V $1.lectins.vcf -sn group1 -sn group2 -sn group3 -sn group4 -sn group5 -sn group6 -sn group19 -sn group24 -o $1.lectins.normals.vcf -env -ef -noTrim
java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T SelectVariants -R ~/genome/genome.fa -V $1.lectins.50kb.vcf -sn group1 -sn group2 -sn group3 -sn group4 -sn group5 -sn group6 -sn group19 -sn group24 -o $1.lectins.normals.50kb.vcf -env -ef -noTrim

#Diseased
java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T SelectVariants -R ~/genome/genome.fa -V $1.lectins.vcf -sn group7 -sn group8 -sn group9 -sn group10 -sn group11 -sn group12 -sn group13 -sn group14 -sn group15 -sn group16 -sn group17 -sn group18 -sn group20 -o $1.lectins.diseased.vcf -ef -env -noTrim
java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T SelectVariants -R ~/genome/genome.fa -V $1.lectins.50kb.vcf -sn group7 -sn group8 -sn group9 -sn group10 -sn group11 -sn group12 -sn group13 -sn group14 -sn group15 -sn group16 -sn group17 -sn group18 -sn group20 -o $1.lectins.diseased.50kb.vcf -ef -env -noTrim

echo .
echo .
echo .

echo done.

echo .
echo .
echo .

#Make the directory structure for the breakdown
echo -n Making the directories for all your pretty new files...
for snpdirs in upstream downstream coding intron all 50kb_up all_by_region
do
if [[ $snpdirs == "all_by_region" ]]
then 
	mkdir $snpdirs
else
	mkdir $snpdirs
	mkdir $snpdirs/mbl1
	mkdir $snpdirs/mbl2
	mkdir $snpdirs/sftpa
	mkdir $snpdirs/sftpd
	mkdir $snpdirs/fcn1
	mkdir $snpdirs/fcn1-like
	mkdir $snpdirs/fcn3
	mkdir $snpdirs/colec10
	mkdir $snpdirs/colec11
	mkdir $snpdirs/colec12
	mkdir $snpdirs/masp1
	mkdir $snpdirs/masp2
fi
done
echo done.


echo .
echo .
echo .



#Start picking out variants.
echo Sorting variants...

#Coding SNPs first
#Coding SNPs don't use the subscript because there are multiple annotations for possible coding SNPs - e.g. nonsense variant, missense variant, etc.
#All coding SNPs for all genes
echo Coding SNPS...

java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'MBL2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'MBL2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'MBL2')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'SFTPA1') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'SFTPA1') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'SFTPA1')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'SFTPD') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'SFTPD') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'SFTPD')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'COLEC10') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'COLEC10') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'COLEC10')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'COLEC11') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'COLEC11') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'COLEC11')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'COLEC12') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'COLEC12') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'COLEC12')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000000436') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000000436') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000000436')) |  ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000024620') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000024620') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000024620')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'FCN3') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'FCN3') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'FCN3')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'MASP1') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'MASP1') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'MASP1')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'MASP2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'MASP2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'MASP2')))" $1.lectins.vcf > coding/$1.all.coding.snps.vcf

java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'MBL2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'MBL2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'MBL2')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'SFTPA1') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'SFTPA1') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'SFTPA1')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'SFTPD') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'SFTPD') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'SFTPD')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'COLEC10') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'COLEC10') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'COLEC10')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'COLEC11') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'COLEC11') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'COLEC11')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'COLEC12') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'COLEC12') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'COLEC12')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000000436') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000000436') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000000436')) |  ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000024620') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000024620') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000024620')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'FCN3') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'FCN3') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'FCN3')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'MASP1') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'MASP1') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'MASP1')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'MASP2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'MASP2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'MASP2')))" $1.lectins.normals.vcf > coding/$1.normals.coding.snps.vcf

java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'MBL2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'MBL2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'MBL2')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$1') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$1') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$1')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'SFTPA1') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'SFTPA1') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'SFTPA1')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'SFTPD') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'SFTPD') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'SFTPD')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'COLEC10') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'COLEC10') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'COLEC10')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'COLEC11') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'COLEC11') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'COLEC11')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'COLEC12') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'COLEC12') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'COLEC12')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000000436') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000000436') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000000436')) |  ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000024620') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000024620') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000024620')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'FCN3') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'FCN3') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'FCN3')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'MASP1') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'MASP1') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'MASP1')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'MASP2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'MASP2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'MASP2')))" $1.lectins.diseased.vcf > coding/$1.diseased.coding.snps.vcf

#MBL1
echo -n MBL1...
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000023001'))" $1.lectins.vcf > coding/mbl1/$1.mbl1.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000023001'))" $1.lectins.normals.vcf > coding/mbl1/$1.mbl1.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000023001'))" $1.lectins.diseased.vcf > coding/mbl1/$1.mbl1.diseased.coding.snps.vcf
echo done.

#MBL2
echo -n MBL2...
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'MBL2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'MBL2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'MBL2'))" $1.lectins.vcf > coding/mbl2/$1.mbl2.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'MBL2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'MBL2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'MBL2'))" $1.lectins.normals.vcf > coding/mbl2/$1.mbl2.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'MBL2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'MBL2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'MBL2'))" $1.lectins.diseased.vcf > coding/mbl2/$1.mbl2.diseased.coding.snps.vcf
echo done.

#SFTPA
echo -n SFTPA...
x=sftpa
X=SFTPA1
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.vcf > coding/$x/$1.$x.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.normals.vcf > coding/$x/$1.$x.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.diseased.vcf > coding/$x/$1.$x.diseased.coding.snps.vcf
echo done.

#SFTPD
echo -n SFTPD...
x=sftpd
X=SFTPD
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.vcf > coding/$x/$1.$x.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.normals.vcf > coding/$x/$1.$x.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.diseased.vcf > coding/$x/$1.$x.diseased.coding.snps.vcf
echo done.

#COLEC10
echo -n COLEC10...
x=colec10
X=COLEC10
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.vcf > coding/$x/$1.$x.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.normals.vcf > coding/$x/$1.$x.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.diseased.vcf > coding/$x/$1.$x.diseased.coding.snps.vcf
echo done.

#COLEC11
echo -n COLEC11...
x=colec11
X=COLEC11
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.vcf > coding/$x/$1.$x.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.normals.vcf > coding/$x/$1.$x.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.diseased.vcf > coding/$x/$1.$x.diseased.coding.snps.vcf
echo done.

#COLEC12
echo -n COLEC12...
x=colec12
X=COLEC12
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.vcf > coding/$x/$1.$x.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.normals.vcf > coding/$x/$1.$x.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.diseased.vcf > coding/$x/$1.$x.diseased.coding.snps.vcf
echo done.

#FCN1
echo -n FCN1...
x=fcn1
X=ENSECAG00000000436
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.vcf > coding/$x/$1.$x.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.normals.vcf > coding/$x/$1.$x.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.diseased.vcf > coding/$x/$1.$x.diseased.coding.snps.vcf
echo done.

#FCN1-like
echo -n FCN1-like...
x=fcn1-like
X=ENSECAG00000024620
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.vcf > coding/$x/$1.$x.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.normals.vcf > coding/$x/$1.$x.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.diseased.vcf > coding/$x/$1.$x.diseased.coding.snps.vcf
echo done.

#FCN3
echo -n FCN3...
x=fcn3
X=FCN3
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.vcf > coding/$x/$1.$x.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.normals.vcf > coding/$x/$1.$x.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.diseased.vcf > coding/$x/$1.$x.diseased.coding.snps.vcf
echo done.

#MASP1
echo -n MASP1...
x=masp1
X=MASP1
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.vcf > coding/$x/$1.$x.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.normals.vcf > coding/$x/$1.$x.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.diseased.vcf > coding/$x/$1.$x.diseased.coding.snps.vcf
echo done.

#MASP2
echo -n MASP2...
x=masp2
X=MASP2
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.vcf > coding/$x/$1.$x.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.normals.vcf > coding/$x/$1.$x.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.lectins.diseased.vcf > coding/$x/$1.$x.diseased.coding.snps.vcf
echo done.

echo Coding SNPs done.


##INTRONS##
#Define variables
dir=intron
type=intron_variant
echo Gathering $dir variants...
. snp_processing_subscript.sh

	#SFTPD - because SFTPD is incorrectly annotated...
	echo -n SFTPD...
	x=sftpd
	X=SFTPD
	java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X') && (POS < 88956555)" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
	java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X') && (POS < 88956555)" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
	java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X') && (POS < 88956555)" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
	echo done.

echo All $dir variants sorted.

##DOWNSTREAM##
#Define Variables
dir=downstream
type=downstream_gene_variant
echo Gathering $dir variants...
. snp_processing_subscript.sh
	#SFTPD - should be fine for downstream.
	echo -n SFTPD...
	x=sftpd
	X=SFTPD
	java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
	java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
	java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
	echo done.
echo All $dir variants sorted.

##UPSTREAM 5kb##
#Define variables
dir=upstream
type=upstream_gene_variant
echo Gathering $dir variants...
. snp_processing_subscript.sh
	#SFTPD - grabbing snps from the 5kb region upstream of the gene.
	echo -n SFTPD...
	x=sftpd
	X=SFTPD
	java -jar ~/java/snpEff/SnpSift.jar filter "(CHROM = 'chr1') && (POS > 88956555) && (POS < 88961555)" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
	java -jar ~/java/snpEff/SnpSift.jar filter "(CHROM = 'chr1') && (POS > 88956555) && (POS < 88961555)" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
	java -jar ~/java/snpEff/SnpSift.jar filter "(CHROM = 'chr1') && (POS > 88956555) && (POS < 88961555)" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
	echo done.


echo All $dir variants sorted.

##UPSTREAM 50kb##
dir=50kb_up
type=upstream_gene_variant
echo Gathering $dir - aka variants within 50kb of the start codon...
. snp_processing_subscript_2.sh
echo All $dir variants sorted.


#ALL BY REGION##
#Grab the total number of SNPs in each region targeted by probes.
dir=all_by_region
echo -n Getting $dir nsps for all genes...
java -jar ~/java/snpEff/SnpSift.jar filter "(((POS > 88892827) && (POS < 89006558)) | ((POS > 43242495) && (POS < 43299033)) | ((POS > 62213524) && (POS < 62301778)) | ((POS >  88547817) && (POS < 88602733)) | ((POS > 40892387) && (POS < 41121801)) | ((POS > 36785097) && (POS < 36875451)) | ((POS > 28409278) && (POS < 28420415)) | ((POS > 25141928) && (POS < 25256814)) | ((POS > 40412424) && (POS < 40446898)))" $1.lectins.vcf > $dir/$1.all.$dir.snps.vcf

java -jar ~/java/snpEff/SnpSift.jar filter "(((POS > 88892827) && (POS < 89006558)) | ((POS > 43242495) && (POS < 43299033)) | ((POS > 62213524) && (POS < 62301778)) | ((POS >  88547817) && (POS < 88602733)) | ((POS > 40892387) && (POS < 41121801)) | ((POS > 36785097) && (POS < 36875451)) | ((POS > 28409278) && (POS < 28420415)) | ((POS > 25141928) && (POS < 25256814)) | ((POS > 40412424) && (POS < 40446898)))" $1.lectins.normals.vcf > $dir/$1.normals.$dir.snps.vcf

java -jar ~/java/snpEff/SnpSift.jar filter "(((POS > 88892827) && (POS < 89006558)) | ((POS > 43242495) && (POS < 43299033)) | ((POS > 62213524) && (POS < 62301778)) | ((POS >  88547817) && (POS < 88602733)) | ((POS > 40892387) && (POS < 41121801)) | ((POS > 36785097) && (POS < 36875451)) | ((POS > 28409278) && (POS < 28420415)) | ((POS > 25141928) && (POS < 25256814)) | ((POS > 40412424) && (POS < 40446898)))" $1.lectins.diseased.vcf > $dir/$1.diseased.$dir.snps.vcf
echo done.

##ALL SNPS FOR EACH GENE##
#SUBSCRIPT NOT USED as the variable "type" is not applicable
#Define Variables
dir=all
#ALL SNPS from the complete VCF file and from normals and diseased subgroups.
echo -n Getting $dir snps for all genes...
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].GENE = 'MBL2') | (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].GENE = 'SFTPA1') | (ANN[*].GENE = 'SFTPD') | (ANN[*].GENE = 'COLEC10') | (ANN[*].GENE = 'COLEC11') | ((ANN[*].GENE has 'COLEC12') | (ANN[*].GENE has 'COLEC12-CETN1')) | (ANN[*].GENE = 'ENSECAG00000000436') | (ANN[*].GENE = 'ENSECAG00000024620') | (ANN[*].GENE = 'FCN3') | (ANN[*].GENE = 'MASP1') | (ANN[*].GENE = 'MASP2'))" $1.lectins.50kb.vcf > $dir/$1.all.$dir.snps.vcf

java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].GENE = 'MBL2') | (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].GENE = 'SFTPA1') |(ANN[*].GENE = 'SFTPD') | (ANN[*].GENE = 'COLEC10') | (ANN[*].GENE = 'COLEC11') | ((ANN[*].GENE has 'COLEC12') | (ANN[*].GENE has 'COLEC12-CETN1')) | (ANN[*].GENE = 'ENSECAG00000000436') | (ANN[*].GENE = 'ENSECAG00000024620') | (ANN[*].GENE = 'FCN3') | (ANN[*].GENE = 'MASP1') | (ANN[*].GENE = 'MASP2'))" $1.lectins.normals.50kb.vcf > $dir/$1.normals.$dir.snps.vcf

java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].GENE = 'MBL2') |  (ANN[*].GENE = 'ENSECAG00000023001') | (ANN[*].GENE = 'SFTPA1') | (ANN[*].GENE = 'SFTPD') | (ANN[*].GENE = 'COLEC10') | (ANN[*].GENE = 'COLEC11') | ((ANN[*].GENE has 'COLEC12') | (ANN[*].GENE has 'COLEC12-CETN1')) | (ANN[*].GENE = 'ENSECAG00000000436') | (ANN[*].GENE = 'ENSECAG00000024620') | (ANN[*].GENE = 'FCN3') | (ANN[*].GENE = 'MASP1') | (ANN[*].GENE = 'MASP2'))" $1.lectins.diseased.50kb.vcf > $dir/$1.diseased.$dir.snps.vcf
echo done.


#MBL1
echo -n MBL1...
x=mbl1
X=ENSECAG00000023001
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#MBL2
echo -n MBL2...
x=mbl2
X=MBL2
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#SFTPA
echo -n SFTPA1...
x=sftpa
X=SFTPA1
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#SFTPD
echo -n SFTPD...
x=sftpd
X=SFTPD
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#COLEC10
echo -n COLEC10...
x=colec10
X=COLEC10
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#COLEC11
echo -n COLEC11...
x=colec11
X=COLEC11
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#COLEC12
echo -n COLEC12...
x=colec12
X=COLEC12
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].GENE has '$X') | (ANN[*].GENE has '$X-'))" $1.lectins.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].GENE has '$X') | (ANN[*].GENE has '$X-'))" $1.lectins.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].GENE has '$X') | (ANN[*].GENE has '$X-'))" $1.lectins.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#FCN1
echo -n FCN1...
x=fcn1
X=ENSECAG00000000436
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#FCN1-like
echo -n FCN1-like...
x=fcn1-like
X=ENSECAG00000024620
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#FCN3
echo -n FCN3...
x=fcn3
X=FCN3
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#MASP1
echo -n MASP1...
x=masp1
X=MASP1
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#MASP2
echo -n MASP2...
x=masp2
X=MASP2
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.lectins.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

echo .
echo .
echo .
echo .
echo .



#Run vcf-compare on all the different files that were just generated
echo -n Generating stats file...
for dir in coding intron downstream upstream all 50kb_up all_by_region
do
	echo Getting stats from $dir variants...
	. snp_stats_subscript.sh

done

#Call in the script that will summarize everything.
snp_stat_summarizer.sh

#Cleanup
mv *50kb.* 50kb_up/

echo DONE.