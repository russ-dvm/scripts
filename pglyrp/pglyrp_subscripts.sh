#!/bin/bash

#RSF

#Used as a subscript in the process_pglyrps.sh script. Requires the variables "dir" and "type" to be defined, should be inherited or defined within this script. 

#gather all types of snps, incl intron, intergenic, upstream, downstream, and 50kb upstream.

for dir in intron upstream downstream 50kb_up
do
	echo $dir
	if [[ $dir != "50kb_up" ]]
	then
		for X in ${!gene_hash[@]}
		do
			# if [[ $X == "ENSECAG00000010847" ]]
			# then
			# 	x=pglyrp1
			# else
			# 	x=$(echo $X | tr [:upper:] [:lower:])
			# fi

			gene_name="${gene_hash[$X]}"


			if [[ $dir == "intron" ]]
			then
				type=intron_variant
			elif [[ $dir == "upstream" ]]
			then
				type=upstream_gene_variant
			else [[ $dir == "downstream" ]]
				type=downstream_gene_variant
			fi
			echo "$gene_name"...
			java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENEID = '$X')" $1.pglyrp.vcf > $dir/"$gene_name"/$1."$gene_name".all.$dir.snps.vcf
			java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENEID = '$X')" $1.pglyrp.normals.vcf > $dir/"$gene_name"/$1."$gene_name".normals.$dir.snps.vcf
			java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENEID = '$X')" $1.pglyrp.diseased.vcf > $dir/"$gene_name"/$1."$gene_name".diseased.$dir.snps.vcf

			#ALL SNPS from the complete VCF file and from normals and diseased subgroups for each gene region.
			echo -n Getting $dir snps for all genes...
			java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYRP2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYR3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYRP4') && (POS < 88956555)))" $1.pglyrp.vcf > $dir/$1.all.$dir.snps.vcf
			java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYRP2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYR3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYRP4') && (POS < 88956555)))" $1.pglyrp.normals.vcf > $dir/$1.normals.$dir.snps.vcf
			java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYRP2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYR3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYRP4') && (POS < 88956555)))" $1.pglyrp.diseased.vcf > $dir/$1.diseased.$dir.snps.vcf
			echo done.	

						
		done
	else
		for X in ${!gene_hash[@]}
		do
			gene_name="${gene_hash[$X]}"

			echo "$gene_name"...
			type=upstream_gene_variant
			java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENEID = '$X')" $1.pglyrp.50kb.vcf > $dir/"$gene_name"/$1."$gene_name".all.$dir.snps.vcf
			java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENEID = '$X')" $1.pglyrp.normals.50kb.vcf > $dir/"$gene_name"/$1."$gene_name".normals.$dir.snps.vcf
			java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENEID = '$X')" $1.pglyrp.diseased.50kb.vcf > $dir/"$gene_name"/$1."$gene_name".diseased.$dir.snps.vcf	

			#ALL SNPS from the complete VCF file and from normals and diseased subgroups for each gene region.
			echo -n Getting $dir snps for all genes...
			java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYRP2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYR3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYRP4') && (POS < 88956555)))" $1.pglyrp.vcf > $dir/$1.all.$dir.snps.vcf
			java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYRP2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYR3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYRP4') && (POS < 88956555)))" $1.pglyrp.normals.vcf > $dir/$1.normals.$dir.snps.vcf
			java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYRP2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYR3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'PGLYRP4') && (POS < 88956555)))" $1.pglyrp.diseased.vcf > $dir/$1.diseased.$dir.snps.vcf
			echo done.	
		done
	fi	
done