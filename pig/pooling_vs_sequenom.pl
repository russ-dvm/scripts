#!/usr/bin/perl

use strict;
use warnings;

#RSF, rfrase03@uoguelph.ca, 2016/06/02
#This script is designed to check the allele frequency of pooled data vs the data obtained via Sequenom genotyping of individual animal within a pool.
#e.g. take a group composed of pig1, pig2, and pig3, with genotype: 0/0/0/1/1/1
#we also have individual data:
# pig1: 0/0
# pig2: 1/1
# pig3: 0/1
# thus, we have 3 ref and 3 alt in both experiments. This script can generate # of ref, # alt for both experiments.
#Data needs a bit of pre-processing - see evernote, but in general, run variantstotable to get letter-form genotypes (A/T instead 0/1), remove extraneous info, and remove slashes using sed.


open (my $vcf, "<", $ARGV[0]);
my $switch = $ARGV[1];

my @vcf_array = <$vcf>;

if ($switch =~ /illumina/){print "group\tid\tref_illumina\talt_illumina\n";}
else {print "group\tid\tref_sequenom\talt_sequenom\n"; }

my (@vcf_header);

foreach my $vcf_line (@vcf_array) {

	chomp $vcf_line;

	if ($vcf_line =~ /ID/) {
		@vcf_header = split(/\t/, $vcf_line);
	}


	else {

		my @vcf_fields = split(/\t/, $vcf_line);
		my $ref_allele = $vcf_fields[1];
		my $alt_allele = $vcf_fields[2];

		for (my $i=3; $i < scalar @vcf_fields; $i++){

			print "$vcf_header[$i]\t$vcf_fields[0]\t";

			if ($vcf_fields[$i] =~ /\.\.\./){
				print "NA\tNA\n";
			}
			else {

				my @ref_count_array = ($vcf_fields[$i] =~ /$ref_allele/g);
				my $ref_count = @ref_count_array;
				my @alt_count_array = ($vcf_fields[$i] =~ /$alt_allele/g);
				my $alt_count = @alt_count_array;
				print "$ref_count\t$alt_count\n";

			}

		}

	}

}