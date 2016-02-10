#!/bin/bash

#RSF 2015/12/03, rfrase03@uoguelph.ca
#Set of commands to select variants from PGLYRPs only; to separate them out by region (upstream, downstream, coding, introns).

#annotate with snpeff
echo -n Running snpEff...
java -Xmx10G -jar ~/java/snpEff/snpEff.jar eff -i vcf -o vcf EquCab2.78 -interval ~/equine/2014_11_24/ref_files/intervals.txt -stats $1.html $1.vcf > $1.snpeff.vcf
java -Xmx10G -jar ~/java/snpEff/snpEff.jar eff -i vcf -o vcf EquCab2.78 -ud 50000 -interval ~/equine/2014_11_24/ref_files/intervals.txt -stats $1.html $1.vcf > $1.snpeff.50kb.vcf
echo done.

#Select variants from pglyrps only -ef = exclude filtered, -env = exclude non variant
echo -n Generating VCFs with variants that are in PGLYRPs and have passed...

java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R ~/genome/genome.fa -L ~/equine/2014_11_24/ref_files/pglyrp_intervals.intervals -ef -env -V $1.snpeff.vcf -o $1.pglyrp.vcf
java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R ~/genome/genome.fa -L ~/equine/2014_11_24/ref_files/pglyrp_intervals.intervals -ef -env -V $1.snpeff.50kb.vcf -o $1.pglyrp.50kb.vcf

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
#5 kb upstream and downstrewam
java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R ~/genome/genome.fa -V $1.pglyrp.vcf -sn group1 -sn group2 -sn group3 -sn group4 -sn group5 -sn group6 -sn group19 -sn group24 -o $1.pglyrp.normals.vcf -env -ef -noTrim

#50kb upstream/downstream
java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R ~/genome/genome.fa -V $1.pglyrp.50kb.vcf -sn group1 -sn group2 -sn group3 -sn group4 -sn group5 -sn group6 -sn group19 -sn group24 -o $1.pglyrp.normals.50kb.vcf -env -ef -noTrim

#Diseased
#5kb upstream/downstream
java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R ~/genome/genome.fa -V $1.pglyrp.vcf -sn group7 -sn group8 -sn group9 -sn group10 -sn group11 -sn group12 -sn group13 -sn group14 -sn group15 -sn group16 -sn group17 -sn group18 -sn group20 -o $1.pglyrp.diseased.vcf -ef -env -noTrim
#50kb upstream/downstream
java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R ~/genome/genome.fa -V $1.pglyrp.50kb.vcf -sn group7 -sn group8 -sn group9 -sn group10 -sn group11 -sn group12 -sn group13 -sn group14 -sn group15 -sn group16 -sn group17 -sn group18 -sn group20 -o $1.pglyrp.diseased.50kb.vcf -ef -env -noTrim

echo .
echo .
echo .

echo done.

echo .
echo .
echo .

#Make the directory structure for the breakdown
echo -n Making the directories for all your pretty new files...
for snpdirs in upstream downstream coding intron all_by_gene 50kb_up all_by_region
do
if [[ $snpdirs == "all_by_region" ]]
then 
	mkdir $snpdirs
else
	mkdir $snpdirs
	mkdir $snpdirs/pglyrp1
	mkdir $snpdirs/pglyrp2
	mkdir $snpdirs/pglyrp3
	mkdir $snpdirs/pglyrp4
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

java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP2')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP3') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP3') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP3')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP4') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP4') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP4')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000010847') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000010847') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000010847')))" $1.pglyrp.vcf > coding/$1.all.coding.snps.vcf

java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP2')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP3') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP3') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP3')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP4') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP4') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP4')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000010847') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000010847') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000010847')))" $1.pglyrp.normals.vcf > coding/$1.normals.coding.snps.vcf

java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP2')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP3') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP3') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP3')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP4') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP4') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP4')) | ((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'ENSECAG00000010847') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'ENSECAG00000010847') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'ENSECAG00000010847')))" $1.pglyrp.diseased.vcf > coding/$1.diseased.coding.snps.vcf

#PGLYRP3
echo -n PGLYRP3...
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP3') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP3') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP3'))" $1.pglyrp.vcf > coding/pglyrp3/$1.pglyrp3.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP3') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP3') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP3'))" $1.pglyrp.normals.vcf > coding/pglyrp3/$1.pglyrp3.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP3') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP3') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP3'))" $1.pglyrp.diseased.vcf > coding/pglyrp3/$1.pglyrp3.diseased.coding.snps.vcf
echo done.

#PGLYRP2
echo -n PGLYRP2...
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP2'))" $1.pglyrp.vcf > coding/pglyrp2/$1.pglyrp2.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP2'))" $1.pglyrp.normals.vcf > coding/pglyrp2/$1.pglyrp2.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = 'PGLYRP2') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = 'PGLYRP2') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = 'PGLYRP2'))" $1.pglyrp.diseased.vcf > coding/pglyrp2/$1.pglyrp2.diseased.coding.snps.vcf
echo done.

#PGLYRP4
echo -n PGLYRP4...
x=pglyrp4
X=PGLYRP4
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.pglyrp.vcf > coding/$x/$1.$x.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.pglyrp.normals.vcf > coding/$x/$1.$x.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.pglyrp.diseased.vcf > coding/$x/$1.$x.diseased.coding.snps.vcf
echo done.

#PGLYRP1
echo -n PGLYRP1...
x=pglyrp1
X=ENSECAG00000010847
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.pglyrp.vcf > coding/$x/$1.$x.all.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.pglyrp.normals.vcf > coding/$x/$1.$x.normals.coding.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'nonsense_variant') && (ANN[*].GENE = '$X') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$X'))" $1.pglyrp.diseased.vcf > coding/$x/$1.$x.diseased.coding.snps.vcf
echo done.


echo Coding SNPs done.


. pglyrp_subscripts.sh


#ALL BY REGION##
#Grab the total number of SNPs in each region targeted by probes.
dir=all_by_region
echo -n Getting $dir nsps for all genes...
java -jar ~/java/snpEff/SnpSift.jar filter "(((POS > 16542802) && (POS < 16585723)) | ((POS > 16445016) && (POS < 16513762)) | ((POS > 704631) && (POS < 728959)) | ((POS >  44235788) && (POS < 44297640)))" $1.pglyrp.vcf > $dir/$1.all.$dir.snps.vcf

java -jar ~/java/snpEff/SnpSift.jar filter "(((POS > 16542802) && (POS < 16585723)) | ((POS > 16445016) && (POS < 16513762)) | ((POS > 704631) && (POS < 728959)) | ((POS >  44235788) && (POS < 44297640)))" $1.pglyrp.normals.vcf > $dir/$1.normals.$dir.snps.vcf

java -jar ~/java/snpEff/SnpSift.jar filter "(((POS > 16542802) && (POS < 16585723)) | ((POS > 16445016) && (POS < 16513762)) | ((POS > 704631) && (POS < 728959)) | ((POS >  44235788) && (POS < 44297640)))" $1.pglyrp.diseased.vcf > $dir/$1.diseased.$dir.snps.vcf
echo done.

##ALL SNPS FOR EACH GENE##
#SUBSCRIPT NOT USED as the variable "type" is not applicable
#Define Variables
dir=all_by_gene
#ALL SNPS from the complete VCF file and from normals and diseased subgroups.
echo -n Getting $dir snps for all genes...
java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].GENE = 'PGLYRP2') | (ANN[*].GENE = 'PGLYRP3') | (ANN[*].GENE = 'PGLYRP4') | (ANN[*].GENE = 'ENSECAG00000010847'))" $1.pglyrp.50kb.vcf > $dir/$1.all.$dir.snps.vcf

java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].GENE = 'PGLYRP2') | (ANN[*].GENE = 'PGLYRP3') | (ANN[*].GENE = 'PGLYRP4') |(ANN[*].GENE = 'ENSECAG00000010847'))" $1.pglyrp.normals.50kb.vcf > $dir/$1.normals.$dir.snps.vcf

java -jar ~/java/snpEff/SnpSift.jar filter "((ANN[*].GENE = 'PGLYRP2') |  (ANN[*].GENE = 'PGLYRP3') | (ANN[*].GENE = 'PGLYRP4') | (ANN[*].GENE = 'ENSECAG00000010847'))" $1.pglyrp.diseased.50kb.vcf > $dir/$1.diseased.$dir.snps.vcf
echo done.


#PGLYRP3
echo -n PGLYRP3...
x=pglyrp3
X=PGLYRP3
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.pglyrp.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.pglyrp.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.pglyrp.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#PGLYRP2
echo -n PGLYRP2...
x=pglyrp2
X=PGLYRP2
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.pglyrp.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.pglyrp.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.pglyrp.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#PGLYRP4
echo -n PGLYRP4...
x=pglyrp4
X=PGLYRP4
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.pglyrp.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.pglyrp.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.pglyrp.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#ENSECAG00000010847
echo -n PGLYRP1...
x=pglyrp1
X=ENSECAG00000010847
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.pglyrp.50kb.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.pglyrp.normals.50kb.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].GENE = '$X')" $1.pglyrp.diseased.50kb.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.


echo .
echo .
echo .
echo .
echo .



#Run vcf-compare on all the different files that were just generated
echo Generating vcf-compare files...
for dir in coding intron downstream upstream all_by_gene 50kb_up all_by_region
do
	echo Comparing $dir variants...
	. pglyrp_stats_subscript.sh

done

#Call in the script that will summarize everything.
pglyrp_stat_summarizer.sh

#Cleanup
mv *50kb.* 50kb_up/

echo DONE.