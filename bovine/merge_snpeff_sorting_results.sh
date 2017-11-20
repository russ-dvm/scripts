#!/bin/bash

## Merge the output of the snpeff_sorting script, since I didn't write in a snpsift command to do it... Weird that GATK doesn't accept lists of VCF files for this function.

java -jar ~/java/GATK/GenomeAnalysisTK.jar \
    -R ~/bovine/genome/genome.chr.fa \
    -T CombineVariants \
    --assumeIdenticalSamples \
    -V cgn/*all*vcf \
    -V colec10/*all*vcf \
    -V colec11/*all*vcf \
    -V colec12/*all*vcf \
    -V colec43/*all*vcf \
    -V colec46/*all*vcf \
    -V fcn1/*all*vcf \
    -V masp1/*all*vcf \
    -V masp2/*all*vcf \
    -V mbl1/*all*vcf \
    -V mbl2/*all*vcf \
    -V sftpa/*all*vcf \
    -V sftpd/*all*vcf \
    -o $1