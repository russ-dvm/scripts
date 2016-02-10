#!/bin/bash

#RSF, rfrase03@uoguelph.ca
#Written to flow from snpeff_sorting.sh


read -p "lectin or pglyrp? if you don't type one things are gonna get real ugly... " family 

if [[ $family == "lectin" ]]
then
	gene_array=( colec10 colec11 colec12 colec43 colec46 cgn mbl1 mbl2 sftpa sftpd fcn1 fcn3 masp1 masp2 )
else
	gene_array=( pglyrp1 pglyrp2 pglyrp3 pglyrp4 )
fi

###FOR WHATEVER REASON OMITS PGLYRP1 THE FIRST TIME THE SCRIPT IS RUN... WTF.

touch $family/summary_stats.txt
echo -e "GENE"'\t'"TYPE"'\t'"TOTAL"'\t'"TOTAL NORMAL"'\t'"UNIQUE NORMAL"'\t'"TOTAL DISEASED"'\t'"UNIQUE DISEASED" >> $family/summary_stats.txt



for type in coding intron upstream downstream 50kb_up
do
	touch $family/$type/"$type"_summary.txt
	echo -e "GENE"'\t'"TYPE"'\t'"TOTAL"'\t'"TOTAL NORMAL"'\t'"UNIQUE NORMAL"'\t'"TOTAL DISEASED"'\t'"UNIQUE DISEASED" >> $family/$type/"$type"_summary.txt	
	for gene in ${gene_array[*]}
	do
		echo summarizing:$type,$gene
# 		all_count=$(grep -vc \# $family/$type/$gene/*$gene.all.$type*)
		norm_count=$(grep -vc \# $family/$type/$gene/*norm*)
		dis_count=$(grep -vc \# $family/$type/$gene/*dis*)
		caps_gene=$(echo $gene | tr [:lower:] [:upper:])
		unique_norm=$(grep ^VN $family/$type/$gene/*compared* | grep normal | grep -v diseased | cut -f 2)
		unique_dis=$(grep ^VN $family/$type/$gene/*compared* | grep diseased | grep -v normal | cut -f 2)
			if [[ -z "$unique_norm" ]] #set empty variables to 0
			then
				unique_norm=0
			fi
			if [[ -z "$unique_dis" ]] 
			then
				unique_dis=0
			fi
		echo -e "$caps_gene"'\t'"$type"'\t'"$all_count"'\t'"$norm_count"'\t'"$unique_norm"'\t'"$dis_count"'\t'"$unique_dis" >> $family/$type/$type_summary.txt
# 		total_gene_count=$(grep -vc \# $family/$type/*all.$type.snp*)
# 		total_norm_count=$(grep -vc \# $family/$type/*normals.$type.snp*)
# 		total_dis_count=$(grep -vc \# $family/$type/*diseased.$type.snp*)
# 		total_unique_norm=$(grep ^VN $family/$type/*compared* | grep normal | grep -v diseased | cut -f 2)
# 		total_unique_dis=$(grep ^VN $family/$type/*compared* | grep diseased | grep -v normal | cut -f 2)
	done
tail -n+2 $family/$type/$type_summary.txt >> $family/summary_stats.txt
done
# echo -e "TOTALS"'\t'"all_by_region"'\t'"$total_gene_count"'\t'"$total_norm_count"'\t'"$total_unique_norm"'\t'"$total_dis_count"'\t'"$total_unique_dis" >> $family/summary_stats.txt	