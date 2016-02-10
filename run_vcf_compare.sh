#simple script to compress vcf files, run vcf compare, and then decompress and remove cluttering index files.
#usage: run_vcf_compare <input file 1> <input file 2> <output>


echo -n Compressing and indexing...
bgzip $1
tabix -p vcf $1.gz
bgzip $2
tabix -p vcf $2.gz
echo done.

echo -n Running VCF Compare...
vcf-compare $1.gz $2.gz > $3
echo done.

echo -n Cleaning up...
rm *tbi
bgzip -d $1.gz
bgzip -d $2.gz
echo done.

