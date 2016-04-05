#!/bin/bash



#@RG Info
#example from library1 fastq file: @D00430:183:C908BANXX:1:1107:1274:2037 1:N:0:CGCATAC
#example from library2 fastq file: @D00430:183:C908BANXX:2:2111:19375:1999 1:N:0:ACAGCAG
#								instrumentID:run#:flowcell_ID:lane:tile:x-pos:y-pos	read:filtered:bit:index

declare -A libhash
source /media/russ/data/porcine/ref_files/libhash.txt

genome=/media/russ/data/porcine/genome/Sequence/Custom/genome.custom.fa
trimmed_dir=/media/russ/data/porcine/2_trimmed
aligned_dir=/media/russ/data/porcine/3_aligned


#use pig or group to process in batch, or remove
for var in pig*R1*fastq.gz
do

# 	lane=$(head -1 $var | cut -d : -f 4)
	base="${var%_R*}"
	LIB="${libhash[$base]}"
	
	if [[ $LIB == "lib1" ]]
	then
		LANE=1
	else
		LANE=2
	fi
	
	read1="$base"_R1_trimmed.fastq.gz
	read2="$base"_R2_trimmed.fastq.gz
	
	ID=C908BANXX."$base".lane"$LANE"
	PU=C908BANXX."$base".lane"$LANE"
	PL=illumina
# 	echo $ID
# 	echo $LIB
# 	echo $PL
	
	SM="$base"
	
	
	bwa mem -k 15 -t 8 -M -R "@RG\tID:"$ID"\tSM:"$SM"\tPL:"$PL"\tPU:"$PU"\tLB:"$LIB"" $genome "$trimmed_dir"/"$read1" "$trimmed_dir"/"$read2" > "$aligned_dir"/"$base"_unsorted.sam
	
	
done

