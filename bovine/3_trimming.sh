#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#trim all fastq files with the same parameters.

uncompressed_dir=/media/russ/data/bovine/2015_12_13/2_uncompressed
trimmed_dir=/media/russ/data/bovine/2015_12_13/3_trimmed
for file in $uncompressed_dir/*r1.fastq
do
	read1="${file#/media/russ/data/bovine/2015_12_13/2_uncompressed/}"
	group_num="${read1%1.fastq}"
	read2="$group_num"2.fastq
	just_num="${group_num%r}"
	java -Xms4g -Xmx4g -jar ~/java/trimmomatic/trimmomatic.jar PE -phred33 "$uncompressed_dir"/"$read1" "$uncompressed_dir"/"$read2" "$trimmed_dir"/"$just_num"r1_trimmed.fastq "$trimmed_dir"/"$just_num"r1_unpaired.fastq "$trimmed_dir"/"$just_num"r2_trimmed.fastq "$trimmed_dir"/"$just_num"r2_unpaired.fastq ILLUMINACLIP:/home/russ/java/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:8:true MAXINFO:70:0.7 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:75 CROP:299
done