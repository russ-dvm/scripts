#!/bin/bash

#RSF

#Used as a subscript in snpeff_sorting.sh

#$base
#

read -p "Refresh my memory - lectin or pglyrp? (TYPE ONE OF THOSE EXACTLY, OR EVERYTHING GOES TO HELL)? " family


for x in coding intron upstream downstream 50kb_up
do
	if [[ $family == "lectin" ]]
	then
		for gene_array in colec10 colec11 colec12 colec43 colec46 cgn mbl1 mbl2 sftpa sftpd fcn1 fcn3 masp1 masp2
		do
			for nd in normal diseased
			do
				bgzip $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf
				tabix -p vcf $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf.gz
			done

			vcf-compare $family/$x/$gene_array/$base.normal.$x.$gene_array.vcf.gz $family/$x/$gene_array/$base.diseased.$x.$gene_array.vcf.gz > $family/$x/$gene_array/$base.$x.$gene_array.compared.txt


	# 		bgzip $dir/$1.diseased.$dir.snps.vcf
	# 		bgzip $dir/$1.normals.$dir.snps.vcf
	# 		tabix -p vcf $dir/$1.diseased.$dir.snps.vcf.gz
	# 		tabix -p vcf $dir/$1.normals.$dir.snps.vcf.gz
	# 		vcf-compare $dir/$1.normals.$dir.snps.vcf.gz $dir/$1.diseased.$dir.snps.vcf.gz > $dir/$1.$dir.snps.compared.txt

		
			for nd in normal diseased
			do
				bgzip -d $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf.gz
			done
		rm $family/$x/$gene_array/*tbi
		done

	else #PGLYRPS
	
		for gene_array in pglyrp1 pglyrp2 pglyrp3 pglyrp4
		do
			for nd in normal diseased
			do
				bgzip $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf
				tabix -p vcf $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf.gz
			done

			vcf-compare $family/$x/$gene_array/$base.normal.$x.$gene_array.vcf.gz $family/$x/$gene_array/$base.diseased.$x.$gene_array.vcf.gz > $family/$x/$gene_array/$base.$x.$gene_array.compared.txt
		
			for nd in normal diseased
			do
				bgzip -d $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf.gz
			done
		rm $family/$x/$gene_array/*tbi
		done
		
	fi
done


