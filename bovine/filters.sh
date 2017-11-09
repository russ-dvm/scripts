#!/usr/bash

# Runs variant filtration on the VCF call files. Separates into INDELS and SNPs and then recombines

genome=/media/russ/data/bovine/genome/Sequence/WholeGenomeFasta/genome.chr.fa

vcf=$1

java -jar ~/java/GATK/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R "$genome" \
    -V $vcf \
    -selectType SNP \
    -o snps.vcf \
    --restrictAllelesTo BIALLELIC 

java -jar ~/java/GATK/GenomeAnalysisTK.jar \
    -R "$genome" \
    -T SelectVariants \
    -V $vcf \
    -selectType INDEL \
    --restrictAllelesTo BIALLELIC \
    -o indel.vcf 


#SNPS
java -jar ~/java/GATK/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R "$genome" \
    -V snps.vcf \
    -filter "QD < 1.5" --filterName "QDFilter" \
    -filter "DP > 890 ? SOR > 3.0 : FS > 60.0" --filterName "SOR-FSFilter" \
    -filter "MQ < 40" --filterName "MQFilter" \
    -G_filter "DP < 50" --genotypeFilterName "DP50" \
    --setFilteredGtToNocall \
    -o snps.filtered.depth50.vcf 

#INDELS
java -Xmx28g -jar ~/java/GATK/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R "$genome" \
    -V indel.vcf \
    -filter "QD < 2.0" \
    --filterName "QDFilter" \
    -filter "FS > 200" \
    --filterName "FSFilter" \
    -filter "ReadPosRankSum > -20.0" \
    --filterName "ReadPosRankSumFilter" \
    -G_filter "DP < 50" \
    --genotypeFilterName "DP50" \
    --setFilteredGtToNocall \
    -o indel.filtered.depth50.vcf


#MERGE SNPS/INDELS back together again
java -Xmx28g -jar ~/java/GATK/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R "$genome" \
    --variant snps.filtered.depth50.vcf \
    --variant indel.filtered.depth50.vcf \
    -assumeIdenticalSamples \
    -o all.unfiltered.variants.vcf

## EMIT ONLY biallelic, variant, and unfiltered snps
## Remove spanning deletions
java -Xmx28g -jar ~/java/GATK/GenomeAnalysisTK.jar \
	-R "$genome" \
	-T SelectVariants \
	-env \
	-ef \
	-V all.unfiltered.variants.vcf \
	-o all.filtered.variants.vcf


java -jar ~/java/snpEff/SnpSift.jar filter -n "(ALT = '*')" all.filtered.variants.vcf > all.filtered.no-span.variants.vcf
