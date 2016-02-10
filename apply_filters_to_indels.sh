#!/bin/bash

#Will run variant filtraiton the specified VCF file (should contain only indels), then export the data to table format (csv), transfer the table to Russ' desktop computer, and then display the stats of the vcf file using the filter_stats.sh script.

#My filters
java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T VariantFiltration -R ~/genome/genome.fa -V indel.vcf -filter "QD < 2.0" --filterName "QDFilter" -filter "SOR > 10.0" --filterName "SORFilter" -filter "FS > 200.0" --filterName "FSFilter" -filter "ReadPosRankSum < -20.0" --filterName "ReadPosFilter" -filter "DP/AC < 10.00" --filterName "DP-AC-Filter" -o $1.vcf 

#Jiao et al filters
#java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T VariantFiltration -R ~/genome/genome.fa -V snps.vcf -filter "DP < 445" --filterName "Depth" -filter "MQ < 20.0" --filterName "MQFilter" -o $1.vcf
#
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
