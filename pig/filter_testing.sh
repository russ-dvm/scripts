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
# 			-filter "ReadPosRankSum < "$rprs"" \
# 			--filterName "ReadPosRankSumFilter" \
# 			-G_filter "DP < 10" \
# 			--genotypeFilterName "DP10" \
# 			--setFilteredGtToNocall \
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
		for gq in 20 30
		do
			java -Xmx28g -jar ~/java/GATK/GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R "$genome" \
			-V $1 \
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
			-G_filter "GQ < "$gq"" \
			--genotypeFilterName "GQ30" \
			--setFilteredGtToNocall \
			-o qd"$qd"_mq"$mq"_gq"$gq".vcf 
		done
	done
done



