#!/bin/bash


trimmed=/work/rfraser/pig/trimmed


for file in *R1*gz
do
	read1="${file%.fastq.gz}"
	group_num="${read1%1}"
	read2="$group_num"2
	just_num="${group_num%_R}"
	java -Xmx60G -jar /work/rfraser/programs/trimmomatic-0.35.jar PE $file "$read2".fastq.gz "$trimmed"/"$read1"_trimmed.fastq.gz "$trimmed"/"$read1"_unpaired.fastq.gz "$trimmed"/"$read2"_trimmed.fastq.gz "$trimmed"/"$read2"_unpaired.fastq.gz MAXINFO:70:0.7 SLIDINGWINDOW:5:20 LEADING:20 TRAILING:20 ILLUMINACLIP:/work/rfraser/programs/adapters/TruSeq3-PE-2.fa:2:30:10:3:true MINLEN:40

done