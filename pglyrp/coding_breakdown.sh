#!/bin/bash

#RSF
#Breakdown coding SNPs by synonymous/nonsynonymous by gene, normal, diseased

touch coding_type_breakdown.txt
echo -e "GENE"'\t'"TYPE"'\t'"TOTAL"'\t'"NORMAL"'\t'"DISEASED" >> coding_type_breakdown.txt

for gene in PGLYRP1a PGLYRP1b PGLYRP1x PGLYRP2 PGLYRP3 PGLYRP4
	do
		# caps_gene=$(echo $gene | tr [:lower:] [:upper:])
		all_syn=$(grep -v \# $gene/*$gene.all.* | grep -ci synonymous)
		all_mis=$(grep -v \# $gene/*$gene.all.* | grep -ci missense)
		norm_syn=$(grep -v \# $gene/*norm*  | grep -ci synonymous)
		norm_mis=$(grep -v \# $gene/*norm* | grep -ci missense)
		dis_syn=$(grep -v \# $gene/*diseased* | grep -ci synonymous)
		dis_mis=$(grep -v \# $gene/*diseased* | grep -ci missense)
			if [[ -z "$all_syn" ]] #set empty variables to 0. there is probably a much better way to do this.
			then
				all_syn=0
			fi
			if [[ -z "$all_mis" ]] 
			then
				all_mis=0
			fi
			if [[ -z "$norm_mis" ]] 
			then
				norm_mis=0
			fi
			if [[ -z "$norm_syn" ]] 
			then
				norm_syn=0
			fi
			if [[ -z "$dis_syn" ]] 
			then
				dis_syn=0
			fi
			if [[ -z "$dis_mis" ]] 
			then
				dis_syn=0
			fi
		echo -e "$gene"'\t'"synonymous"'\t'"$all_syn"'\t'"$norm_syn"'\t'"$dis_syn" >> coding_type_breakdown.txt
		echo -e "$gene"'\t'"missense"'\t'"$all_mis"'\t'"$norm_mis"'\t'"$dis_mis" >> coding_type_breakdown.txt
done
echo done.