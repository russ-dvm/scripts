#!/bin/bash

#run vcf-compare without having to compress & index files first

#vcfcompare.sh <file1> <file2> <output>

bgzip $1
tabix -p vcf $1.gz

bgzip $2
tabix -p vcf $2.gz

vcf-compare $1.gz $2.gz > $3

bgzip -d $1.gz
bgzip -d $2.gz

rm $1.gz.tbi
rm $2.gz.tbi