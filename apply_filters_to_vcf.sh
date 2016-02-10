#!/bin/bash

#Will run variant filtraiton on the specified VCF file, then export the data to table format (csv), transfer the table to Russ' desktop computer, and then display the stats of the vcf file using the filter_stats.sh script.

#My filters
java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T VariantFiltration -R ~/genome/genome.fa -V snps.vcf -filter "QD < 1.35" --filterName "QDFilter" -filter "DP > 890 ?	SOR > 2.0 : FS > 40.0" --filterName "SOR-FSFilter" -filter "MQ < 40.0" --filterName "MQFilter" -filter "DP < 445" --filterName "DepthFilter" -o $1.vcf 

#Jiao et al filters
#java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T VariantFiltration -R ~/genome/genome.fa -V snps.vcf -filter "DP < 445" --filterName "Depth" -filter "MQ < 20.0" --filterName "MQFilter" -o $1.vcf

echo .
echo .
echo .
echo .

java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -R ~/genome/genome.fa -T VariantsToTable -V $1.vcf -F CHROM -F POS -F ID -F QUAL -F FILTER -F AC -F AF -F DP -F QD -F SOR -F FS -F MQ -F ReadPosRankSum -o $1.table --allowMissingData --showFiltered -SMA

echo .
echo .
echo .
echo .


scp $1.table Rick@131.104.118.154:~/Desktop/remote/filterTests

echo .
echo .
echo .
echo .


bash filter_stats.sh $1.vcf
