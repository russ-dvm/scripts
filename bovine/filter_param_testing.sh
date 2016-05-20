#!/bin/bash

#Will run variant filtraiton on the specified VCF file, then export the data to table format (csv), transfer the table to Russ' desktop computer, and then display the stats of the vcf file using the filter_stats.sh script.


genome=/media/russ/data/bovine/genome/Sequence/WholeGenomeFasta/genome.chr.fa


# #My SNP filters
# for qd in 0.5 1.0 1.5 2.0
# do
# 	for mq in 20.0 30.0 40.0
# 	do
# 		java -jar ~/java/GATK/GenomeAnalysisTK.jar -T VariantFiltration -R "$genome" -V snps.vcf -filter "QD < $qd" --filterName "QDFilter" -filter "DP > 890 ? SOR > 2.0 : FS > 40.0" --filterName "SOR-FSFilter" -filter "MQ < $mq" --filterName "MQFilter" -G_filter "DP < 50" --genotypeFilterName "DP50" --setFilteredGtToNocall -o qd"$qd".mq"$mq".depth50.vcf 
# 	done
# done


## My INDEL filters
for qd in 1.0 2.0
do
	for fs in 175.0 200.0 225.0
	do
		for rprs in -15.0 -20.0 -25.0
		do
			java -Xmx28g -jar ~/java/GATK/GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R "$genome" \
			-V $1 \
			-filter "QD < "$qd"" \
			--filterName "QDFilter" \
			-filter "FS > "$fs"" \
			--filterName "FSFilter" \
			-filter "ReadPosRankSum > "$rprs"" \
			--filterName "ReadPosRankSumFilter" \
			-G_filter "DP < 50" \
			--genotypeFilterName "DP50" \
			--setFilteredGtToNocall \
			-o indel.qd"$qd".fs"$fs".vcf
		done
	done
done

#Jiao et al filters
#java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T VariantFiltration -R ~/genome/genome.fa -V snps.vcf -filter "DP < 445" --filterName "Depth" -filter "MQ < 20.0" --filterName "MQFilter" -o $1.vcf

# echo .
# echo .
# echo .
# echo .

# java -jar ~/java/GATK/GenomeAnalysisTK.jar -R "$genome" -T VariantsToTable -V $1.vcf -F CHROM -F POS -F ID -F QUAL -F FILTER -F AC -F AF -F DP -F QD -F SOR -F FS -F MQ -F ReadPosRankSum -o $1.table --allowMissingData --showFiltered -SMA

# echo .
# echo .
# echo .
# echo .


# # scp $1.table Rick@131.104.118.154:~/Desktop/remote/filterTests

# echo .
# echo .
# echo .
# echo .


# bash filter_stats.sh $1.vcf
