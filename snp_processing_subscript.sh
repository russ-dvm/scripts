#!/bin/bash

#RSF

#Used as a subscript in the processing_snps_by_region.sh script. Requires the variables "dir" and "type" to be defined. 

#MBL1
echo -n MBL1...
x=mbl1
X=ENSECAG00000023001
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#MBL2
echo -n MBL2...
x=mbl2
X=MBL2
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#SFTPA
echo -n SFTPA1...
x=sftpa
X=SFTPA1
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#SFTPD
#SFTPD IS ABSENT BECAUSE IT IS BADLY ANNOTATED IN ENSEMBL

#COLEC10
echo -n COLEC10...
x=colec10
X=COLEC10
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#COLEC11
echo -n COLEC11...
x=colec11
X=COLEC11
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#COLEC12
echo -n COLEC12...
x=colec12
X=COLEC12
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#FCN1
echo -n FCN1...
x=fcn1
X=ENSECAG00000000436
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#FCN1-like
echo -n FCN1-like...
x=fcn1-like
X=ENSECAG00000024620
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#FCN3
echo -n FCN3...
x=fcn3
X=FCN3
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#MASP1
echo -n MASP1...
x=masp1
X=MASP1
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#MASP2
echo -n MASP2...
x=masp2
X=MASP2
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.vcf > $dir/$x/$1.$x.all.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.normals.vcf > $dir/$x/$1.$x.normals.$dir.snps.vcf
java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$X')" $1.lectins.diseased.vcf > $dir/$x/$1.$x.diseased.$dir.snps.vcf
echo done.

#ALL SNPS from the complete VCF file and from normals and diseased subgroups.
echo -n Getting $dir snps for all genes...

#DO some fancy juggling because SPD is stupid and badly annotated in Ensembl.
if [[ $type == "intron_variant" ]]
then
	java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MBL2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPA1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPD') && (POS < 88956555)) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC10')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC11')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC12')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000000436')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000024620')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'FCN3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP2')))" $1.lectins.vcf > $dir/$1.all.$dir.snps.vcf
	java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MBL2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPA1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPD') && (POS < 88956555)) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC10')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC11')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC12')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000000436')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000024620')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'FCN3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP2')))" $1.lectins.normals.vcf > $dir/$1.normals.$dir.snps.vcf
	java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MBL2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPA1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPD') && (POS < 88956555)) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC10')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC11')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC12')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000000436')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000024620')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'FCN3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP2')))" $1.lectins.diseased.vcf > $dir/$1.diseased.$dir.snps.vcf
	echo done.
elif [[ $type == "upstream_gene_variant" ]]
then
	java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MBL2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPA1')) | ((CHROM = 'chr1') && (POS > 88956555) && (POS < 88961555)) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC10')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC11')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC12')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000000436')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000024620')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'FCN3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP2')))" $1.lectins.vcf > $dir/$1.all.$dir.snps.vcf
	java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MBL2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPA1')) | ((CHROM = 'chr1') && (POS > 88956555) && (POS < 88961555)) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC10')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC11')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC12')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000000436')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000024620')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'FCN3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP2')))" $1.lectins.normals.vcf > $dir/$1.normals.$dir.snps.vcf
	java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MBL2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPA1')) | ((CHROM = 'chr1') && (POS > 88956555) && (POS < 88961555)) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC10')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC11')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC12')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000000436')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000024620')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'FCN3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP2')))" $1.lectins.diseased.vcf > $dir/$1.diseased.$dir.snps.vcf
	echo done.
else
	java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MBL2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPA1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPD')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC10')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC11')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC12')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000000436')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000024620')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'FCN3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP2')))" $1.lectins.vcf > $dir/$1.all.$dir.snps.vcf
	java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MBL2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPA1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPD')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC10')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC11')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC12')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000000436')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000024620')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'FCN3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP2')))" $1.lectins.normals.vcf > $dir/$1.normals.$dir.snps.vcf
	java -jar ~/java/snpEff/SnpSift.jar filter "(((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MBL2')) |  ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000023001')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPA1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'SFTPD')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC10')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC11')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'COLEC12')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000000436')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'ENSECAG00000024620')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'FCN3')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP1')) | ((ANN[*].EFFECT has '$type') && (ANN[*].GENE = 'MASP2')))" $1.lectins.diseased.vcf > $dir/$1.diseased.$dir.snps.vcf
	echo done.
fi
