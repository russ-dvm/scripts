#!/usr/bin/perl

use strict;
use warnings;

#RFraser, rfrase03@uoguelph.ca
#Quick script to count the number of genotyped versus NA (./.) in a VCF file. Redirect stout to file of your choice.

open (my $vcf, "<", $ARGV[0]);

die "\nPlease specify a VCF file.\nUsage: perl genotype_count.pl <vcf_file> > <output>\n\n" if @ARGV != 1;

my @vcf_file = <$vcf>;

print "location\tnum.genotyped\tnum.na\n";


foreach my $vcf_line (@vcf_file){
	chomp $vcf_line;

	my $genotyped_counter = 0;
	my $na_counter = 0;


	if ($vcf_line =~ /#/){
	}
	
	else {
		my @vcf_fields = split("\t", $vcf_line);
		print "$vcf_fields[0].$vcf_fields[1]\t";

		for (my $i=9; $i <= 77; $i++) {
			if ($vcf_fields[$i] =~ /\.\/\./) {
				$na_counter = $na_counter + 1;
			}
			else {
				$genotyped_counter = $genotyped_counter + 1;

			}

		}

	print "$genotyped_counter\t$na_counter\n";

	}

}