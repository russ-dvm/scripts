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
	"." => "0",
	);
	
die "\nProper usage: perl vcf-to-linkage.pl <vcf-file>. Redirect to stout file of your choice.\n\n" if @ARGV != 1;	
	

foreach my $vcf_line (@vcf_array) {
	
	chomp $vcf_line;

	my @vcf_fields = split(/\t/, $vcf_line);

	if ($vcf_line =~ /##/) {}

	#create header lines by piggy-backing of off the single occurance of #CHROM
	elsif ($vcf_line =~ /#CHROM/) {

		#family info - should be unique for every individual in our case
		print "family\t";
		for (my $i = 0; $i <= $#vcf_fields; $i++) {
			print "$i\t";
		}
		print "\n";

		#individual info - use the original #CHROM header line to obtain pigID
		print "individual_id\t";
		for (my $i = 9; $i <= $#vcf_fields; $i++) {
			print "$vcf_fields[$i]\t";
		}
		print "\n";

		#father_id - unknown, print 0 for all
		print "father_id\t";
		for (my $i = 0; $i <= $#vcf_fields; $i++) {
			print "0\t";
		}
		print "\n";

		#mother_id - unknown, print 0 for all
		print "mother_id\t";
		for (my $i = 0; $i <= $#vcf_fields; $i++) {
			print "0\t";
		}
		print "\n";

		#gender - sort of known. Does it matter? Check with BL.
		print "gender\t";
		for (my $i = 0; $i <= $#vcf_fields; $i++) {
			print "1\t";
		}
		print "\n";

		#disease status - sort of known. Does it matter? Check with BL.
		print "status\t";
		for (my $i = 0; $i <= $#vcf_fields; $i++) {
			print "1\t";
		}
		print "\n";
	}

	else {
	
		
		my $ref_allele = $vcf_fields[3];
		my $alt_allele = $vcf_fields[4];
		
		print "$vcf_fields[0].$vcf_fields[1].$vcf_fields[2]\t";

		for (my $i = 9; $i <= $#vcf_fields; $i++){
		
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
		