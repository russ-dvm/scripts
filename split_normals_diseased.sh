#!/bin/bash
#RSF 2015/10/06
#Quickly split a VCF file containing all groups into two files, one for normals and one for diseased.
#Flags used: -ef = exclude filtered; -env = exclude non variants; noTrim = don't trim alleles (not sure what that means)

#Normals
java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T SelectVariants -R ~/genome/genome.fa -V $1 -sn group1 -sn group2 -sn group3 -sn group4 -sn group5 -sn group6 -sn group19 -sn group24 -o normals.$1 -env -ef -noTrim

#Diseased
java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T SelectVariants -R ~/genome/genome.fa -V $1 -sn group7 -sn group8 -sn group9 -sn group10 -sn group11 -sn group12 -sn group13 -sn group14 -sn group15 -sn group16 -sn group17 -sn group18 -sn group20 -o diseased.$1 -ef -env -noTrim

