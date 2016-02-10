#!/bin/bash

#R.S.F.
#Run the various commands required to generate an allelic association test performed by Plink2.0


echo -n Making PLINK files...
plink2 --vcf $1.with.ids.vcf --keep-allele-order --make-bed --allow-no-sex --make-pheno phenotype.txt 5 --out $1.plink 
echo done.

#echo -n Assinging phenotypes to PLINK files...
#plink2 --bfile $1.plink --make-bed --keep-allele-order --allow-no-sex --a2-allele $1.plink.bim~ --make-pheno phenotype.txt 5  --out $1.plink

echo -n Running allelic association test...
plink2 --bfile $1.plink --model --keep-allele-order --a2-allele $1.plink.bim --allow-no-sex --out $1.model
echo done.

echo -n Cleaning up the MODEL file...
head -n 1 $1.model.model > $1.temp
awk '$5 == "ALLELIC"' $1.model.model >> $1.temp
cat $1.temp | tr -s ' ' '\t' > $1.alleles.txt
rm $1.temp
echo done.

echo -n Sending to remote machine...
scp $1.alleles.txt Rick@131.104.118.154:~/Desktop/remote
echo done.