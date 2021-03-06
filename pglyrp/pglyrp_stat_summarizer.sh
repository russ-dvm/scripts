#!/bin/bash

#RSF
#Generate a TXT file with a summary of all the SNP stats.

touch summary_stats.txt
echo -e "GENE"'\t'"TYPE"'\t'"TOTAL"'\t'"TOTAL_NORMAL"'\t'"UNIQUE_NORMAL"'\t'"TOTAL_DISEASED"'\t'"UNIQUE_DISEASED" >> summary_stats.txt

for type in all_by_gene coding intron upstream downstream 50kb_up
do
	touch $type/$type.txt
	echo -e "GENE"'\t'"TYPE"'\t'"TOTAL"'\t'"TOTAL_NORMAL"'\t'"UNIQUE_NORMAL"'\t'"TOTAL_DISEASED"'\t'"UNIQUE_DISEASED" >> $type/$type.txt	
	for X in ${!gene_hash[@]}
	do
		gene="${gene_hash[$X]}"
		all_count=$(grep -vc \# $type/$gene/*$gene.all.$type*)
		norm_count=$(grep -vc \# $type/$gene/*norm*)
		dis_count=$(grep -vc \# $type/$gene/*dis*)
		# caps_gene=$(echo $gene | tr [:lower:] [:upper:])
		unique_norm=$(grep ^VN $type/$gene/*compared* | grep normal | grep -v diseased | cut -f 2)
		unique_dis=$(grep ^VN $type/$gene/*compared* | grep diseased | grep -v normal | cut -f 2)
			if [[ -z "$unique_norm" ]] #set empty variables to 0
			then
				unique_norm=0
			fi
			if [[ -z "$unique_dis" ]] 
			then
				unique_dis=0
			fi
		echo -e "$gene"'\t'"$type"'\t'"$all_count"'\t'"$norm_count"'\t'"$unique_norm"'\t'"$dis_count"'\t'"$unique_dis" >> $type/$type.txt
		total_gene_count=$(grep -vc \# $type/*all.$type.snp*)
		total_norm_count=$(grep -vc \# $type/*normals.$type.snp*)
		total_dis_count=$(grep -vc \# $type/*diseased.$type.snp*)
		total_unique_norm=$(grep ^VN $type/*compared* | grep normal | grep -v diseased | cut -f 2)
		total_unique_dis=$(grep ^VN $type/*compared* | grep diseased | grep -v normal | cut -f 2)
	done
tail -n+2 $type/$type.txt >> summary_stats.txt
done
# echo -e "TOTALS"'\t'"all_by_region"'\t'"$total_gene_count"'\t'"$total_norm_count"'\t'"$total_unique_norm"'\t'"$total_dis_count"'\t'"$total_unique_dis" >> summary_stats.txt	