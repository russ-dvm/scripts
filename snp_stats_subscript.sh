#!/bin/bash

#RSF

#Written to be used in conjunction with process_snps_by_region.sh . Requires the input of variable "dir".
#Will generate vcf-compare files for all the files created by process_snps_by_region.sh


#Generate vcf-compare files
if [[ $dir == "all_by_region" ]]
then
	bgzip $dir/$1.normals.$dir.snps.vcf
	tabix -p vcf $dir/$1.normals.$dir.snps.vcf.gz
	bgzip $dir/$1.diseased.$dir.snps.vcf
	tabix -p vcf $dir/$1.diseased.$dir.snps.vcf.gz
	vcf-compare $dir/$1.normals.$dir.snps.vcf.gz $dir/$1.diseased.$dir.snps.vcf.gz > $dir/$1.$dir.snps.compared.txt
	rm $dir/*tbi
	bgzip -d $dir/$1.diseased.$dir.snps.vcf.gz
	bgzip -d $dir/$1.normals.$dir.snps.vcf.gz
	echo done.

else
	for x in mbl1 mbl2 sftpa sftpd fcn1 fcn1-like fcn3 colec10 colec11 colec12 masp1 masp2

	do
		echo -n $x...
		bgzip $dir/$x/$1.$x.normals.$dir.snps.vcf
		tabix -p vcf $dir/$x/$1.$x.normals.$dir.snps.vcf.gz
		bgzip $dir/$x/$1.$x.diseased.$dir.snps.vcf
		tabix -p vcf $dir/$x/$1.$x.diseased.$dir.snps.vcf.gz
		vcf-compare $dir/$x/$1.$x.normals.$dir.snps.vcf.gz $dir/$x/$1.$x.diseased.$dir.snps.vcf.gz > $dir/$x/$1.$x.$dir.compared.txt
		bgzip $dir/$1.diseased.$dir.snps.vcf
		bgzip $dir/$1.normals.$dir.snps.vcf
		tabix -p vcf $dir/$1.diseased.$dir.snps.vcf.gz
		tabix -p vcf $dir/$1.normals.$dir.snps.vcf.gz
		vcf-compare $dir/$1.normals.$dir.snps.vcf.gz $dir/$1.diseased.$dir.snps.vcf.gz > $dir/$1.$dir.snps.compared.txt
		rm $dir/$x/*tbi $dir/*tbi
		bgzip -d $dir/$x/$1.$x.normals.$dir.snps.vcf.gz
		bgzip -d $dir/$x/$1.$x.diseased.$dir.snps.vcf.gz
		bgzip -d $dir/$1.diseased.$dir.snps.vcf.gz
		bgzip -d $dir/$1.normals.$dir.snps.vcf.gz
		echo done.
done
fi
