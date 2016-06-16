#!/bin/bash



genome=/work/rfraser/equine/equcab2/Sequence/WholeGenomeFasta/genome.fa
dbsnp=/work/rfraser/equine/snp/dbSNP/merged_chr_dbSNP.vcf
ploidy10=( group1 group2 group24 group10 group13 group14 )
ploidy8=( group3 group4 group5 group6 group19 group7 group8 group9 group11 group12 group20 group15 group16 group17 )
ploidy6=group18
intervals=/work/rfraser/equine/ref_files/intervals.intervals



for x in *.bam
do
	base="${x#/work/rfraser/equine/2014_11_24/realigned/}"
	base="${base%_realigned_recalibrated.bam}"
	
	#set ploidy for different groups
	for y in ${ploidy10[*]}
	do
		if [[ "$y" == "$base" ]]
		then
			ploidy=10
		fi
	done
	for y in ${ploidy8[*]}
	do
		if [[ "$y" == "$base" ]]
		then
		ploidy=8
		fi
	done
	if [[ "$base" == "group18" ]]
	then
		ploidy=6
	fi

# 	echo $base $ploidy



# 	java -Xmx120g -jar /work/rfraser/programs/GenomeAnalysisTK.jar -T HaplotypeCaller -R $genome -ERC GVCF --genotyping_mode DISCOVERY -pcrModel CONSERVATIVE -stand_call_conf 30 -stand_emit_conf 10 -ploidy $ploidy -L $intervals -ip 100 -D $dbsnp -I $x -o ../gvcf/$base.g.vcf -bamout ../gvcf/$base.bamout.bam -maxAltAlleles 4 -log $base.log -maxNumPLValues 20000 -l DEBUG
	
	
	sqsub -r 3.5d -n 8 -mpp 120G -o haplotype.caller.log java -Xmx120g -jar /work/rfraser/programs/GenomeAnalysisTK.jar -T HaplotypeCaller -R $genome -ERC GVCF --genotyping_mode DISCOVERY -pcrModel CONSERVATIVE -stand_call_conf 30 -stand_emit_conf 10 -ploidy $ploidy -L $intervals -ip 100 -D $dbsnp -I $x -o ../gvcf/$base.g.vcf -bamout ../gvcf/$base.bamout.bam -maxAltAlleles 4 -log ../gvcf/$base.log
	

done
