#!/bin/bash

#run the multiple scripts required to set up the VCF files for plink.
#input: plink_setup.sh (filename-omit-extension)

echo -n Filter out failed variants...
java -jar /usr/local/GATK/GenomeAnalysisTK-new.jar -T SelectVariants -R ~/genome/genome.fa -ef -V $1.vcf -o $1.pass.vcf


echo -n "Making all genotypes diploid and removing sample FORMAT information..."
make_diploid.pl $1.pass.vcf $1.diploid.vcf
echo done.

#Cleanup stage no longer needed with modification of the make_diploid script.
#echo -n Cleaning up extraneous information...
#cleanup.pl $1.diploid.vcf $1.cleaned.vcf
#echo done.

echo -n Assigning all variants an ID...
assign.snp.id.pl $1.diploid.vcf $1.with.ids.vcf
echo done.

echo -n "Fixing the header to reflect the appropriate number of diploid animals..."
sed -i 's/#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	group1	group10	group11	group12	group13	group14	group15	group16	group17	group18	group19	group2	group20	group24	group3	group4	group5	group6	group7	group8	group9/#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	group1.group1.1	group1.group1.2	group1.group1.3	group1.group1.4	group1.group1.5	group10.group10.1	group10.group10.2	group10.group10.3	group10.group10.4	group10.group10.5	group11.group11.1	group11.group11.2	group11.group11.3	group11.group11.4	group12.group12.1	group12.group12.2	group12.group12.3	group12.group12.4	group13.group13.1	group13.group13.2	group13.group13.3	group13.group13.4	group13.group13.5	group14.group14.1	group14.group14.2	group14.group14.3	group14.group14.4	group14.group14.5	group15.group15.1	group15.group15.2	group15.group15.3	group15.group15.4	group16.group16.1	group16.group16.2	group16.group16.3	group16.group16.4	group17.group17.1	group17.group17.2	group17.group17.3	group17.group17.4	group18.group18.1	group18.group18.2	group18.group18.3	group19.group19.1	group19.group19.2	group19.group19.3	group19.group19.4	group2.group2.1	group2.group2.2	group2.group2.3	group2.group2.4	group2.group2.5	group20.group20.1	group20.group20.2	group20.group20.3	group20.group20.4	group24.group24.1	group24.group24.2	group24.group24.3	group24.group24.4	group24.group24.5	group3.group3.1	group3.group3.2	group3.group3.3	group3.group3.4	group4.group4.1	group4.group4.2	group4.group4.3	group4.group4.4	group5.group5.1	group5.group5.2	group5.group5.3	group5.group5.4	group6.group6.1	group6.group6.2	group6.group6.3	group6.group6.4	group7.group7.1	group7.group7.2	group7.group7.3	group7.group7.4	group8.group8.1	group8.group8.2	group8.group8.3	group8.group8.4	group9.group9.1	group9.group9.2	group9.group9.3	group9.group9.4/' $1.with.ids.vcf
echo done.