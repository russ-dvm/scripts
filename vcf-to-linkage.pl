#!/usr/bin/perl

use strict;
use warnings;

open (my $vcf, "<", $ARGV[0]);

my @vcf_array = <$vcf>;

my %conversion = (
	A => "1",
	C => "2",
	G => "3",
	T => "4",
	"." => ".",
	);
	
die "\nProper usage: perl vcf-to-linkage.pl <vcf-file>. Redirect to stout file of your choice.\n\n" if @ARGV != 1;	
	
foreach my $vcf_line (@vcf_array) {
	
	chomp $vcf_line;
	
	if ($vcf_line =~ /##/) {}
	elsif ($vcf_line =~ /#CHROM/) {print "$vcf_line\n"}
	else {
	
		my @vcf_fields = split(/\t/, $vcf_line);
		
		my $ref_allele = $vcf_fields[3];
		my $alt_allele = $vcf_fields[4];
		
		print "$vcf_fields[0].$vcf_fields[1].$vcf_fields[2]\t";

		
		for (my $i = 9; $i < $#vcf_fields; $i++){
		
			my @individual_fields = split(/\:/, $vcf_fields[$i]);
			my @individual_genotype = split(/\//, $individual_fields[0]);
			
			my $individual_allele_one = $individual_genotype[0];
			my $individual_allele_two = $individual_genotype[1];
			
			if ($individual_allele_one eq ".")  {
				$individual_allele_one = ".";
			}
			
			elsif ($individual_allele_one == 0) {
				$individual_allele_one = $ref_allele;
			}
			elsif ($individual_allele_one == 1) {
				$individual_allele_one = $alt_allele;
			}

			if ($individual_allele_two eq ".")  {
				$individual_allele_two = ".";
			}

			elsif ($individual_allele_two == 0) {
				$individual_allele_two = $ref_allele;
			}
			elsif ($individual_allele_two == 1) {
				$individual_allele_two = $alt_allele;
			}

			print "$conversion{$individual_allele_one} $conversion{$individual_allele_two}\t";

			}
			print "\n";
		}
	}
		