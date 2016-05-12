#!/bin/bash

#Will run variant filtraiton on the specified VCF file, then export the data to table format (csv), transfer the table to Russ' desktop computer, and then display the stats of the vcf file using the filter_stats.sh script.


genome=/media/russ/data/porcine/genome/genome.custom.fa


#My SNP filters - individuals
# for qd in 1.0 2.0 2.5 3.0
# do
# 	for mq in 30.0 40.0
# 	do
# 		java -Xmx28g -jar ~/java/GATK/GenomeAnalysisTK.jar \
# 		-T VariantFiltration \
# 		-R "$genome" \
# 		-V pig.snps.only.vcf \
# 		-filter "QD < "$qd"" \
# 		--filterName "QDFilter" \
# 		-filter "DP > 720 ? SOR > 2.0 : FS > 40.0" \
# 		--filterName "SOR-FSFilter" \
# 		-filter "MQRankSum < -12.5" \
# 		--filterName "MQRankSumFilter" \
# 		-filter "ReadPosRankSum < -8.0" \
# 		--filterName "ReadPosRankSumFilter" \
# 		-filter "MQ < "$mq"" \
# 		--filterName "MQFilter" \
# 		-G_filter "DP < 10" \
# 		--genotypeFilterName "DP10" \
# 		--setFilteredGtToNocall \
# 		-o qd"$qd"_mq"$mq".vcf 
# 	done
# done


#My INDEL filters
# for qd in 1.0 2.0
# do
# 	for fs in 200.0 175.0 225.0
# 	do
# 		for rprs in -15.0 -20.0 -25.0
# 		do
# 			java -Xmx28g -jar ~/java/GATK/GenomeAnalysisTK.jar \
# 			-T VariantFiltration \
# 			-R "$genome" \
# 			-V pig.indel.only.vcf \
# 			-filter "QD < "$qd"" \
# 			--filterName "QDFilter" \
# 			-filter "FS > "$fs"" \
# 			--filterName "FSFilter" \
# 			-filter "ReadPosRankSum > "$rprs"" \
# 			--filterName "ReadPosRankSumFilter" \
# 			-G_filter "DP < 10" \
# 			--genotypeFilterName "DP10" \
# 		--setFilteredGtToNocall \
# 			-o indel.qd"$qd".fs"$fs".vcf
# 		done
# 	done
# done

###########GROUP FILTERS###########
#My SNP filters - groups
for qd in 1.0 2.0 2.5
do
	for mq in 30.0 40.0
	do
		java -Xmx28g -jar ~/java/GATK/GenomeAnalysisTK.jar \
		-T VariantFiltration \
		-R "$genome" \
		-V group.all.snps.vcf \
		-filter "QD < "$qd"" \
		--filterName "QDFilter" \
		-filter "DP > 720 ? SOR > 2.0 : FS > 40.0" \
		--filterName "SOR-FSFilter" \
		-filter "MQRankSum < -12.5" \
		--filterName "MQRankSumFilter" \
		-filter "ReadPosRankSum < -8.0" \
		--filterName "ReadPosRankSumFilter" \
		-filter "MQ < "$mq"" \
		--filterName "MQFilter" \
		-G_filter "DP < 30" \
		--genotypeFilterName "DP30" \
		--setFilteredGtToNocall \
		-o groups/qd"$qd"_mq"$mq".vcf 
	done
done


#Jiao et al filters
#java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T VariantFiltration -R ~/genome/genome.fa -V snps.vcf -filter "DP < 445" --filterName "Depth" -filter "MQ < 20.0" --filterName "MQFilter" -o $1.vcf

# echo .
# echo .
# echo .
# echo .
# 
# java -jar ~/java/GATK/GenomeAnalysisTK.jar -R "$genome" -T VariantsToTable -V $1.vcf -F CHROM -F POS -F ID -F QUAL -F FILTER -F AC -F AF -F DP -F QD -F SOR -F FS -F MQ -F ReadPosRankSum -o $1.table --allowMissingData --showFiltered -SMA
# 
# echo .
# echo .
# echo .
# echo .
# 
# 
# # scp $1.table Rick@131.104.118.154:~/Desktop/remote/filterTests
# 
# echo .
# echo .
# echo .
# echo .
# 
# 
# bash filter_stats.sh $1.vcf
