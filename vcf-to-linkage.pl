#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Array::Transpose;

open (my $vcf, "<", $ARGV[0]);
my $ped = "$ARGV[1].ped";
my $info = "$ARGV[1].info";
open (my $ped_output, ">", $ped);
open (my $info_output, ">", $info);
my @vcf_array = <$vcf>;

my (@mega_array, @family_array, @individ_array, @father_array, @mother_array, @gender_array, @disease_array);

my %conversion = (
	A => "1",
	C => "2",
	G => "3",
	T => "4",
	"." => "0",
	);
	
die "\nProper usage: perl vcf-to-linkage.pl <vcf-file> <base-name-of-output>. Redirect to stout file of your choice.\n\n" if @ARGV != 2;	
	

foreach my $vcf_line (@vcf_array) {
	
	chomp $vcf_line;

	my @vcf_fields = split(/\t/, $vcf_line);

	if ($vcf_line =~ /##/) {}

	#create header lines by piggy-backing of off the single occurance of #CHROM
	elsif ($vcf_line =~ /#CHROM/) {

		#family info - should be unique for every individual in our case
		# print "family\t";
		for (my $i = 9; $i <= $#vcf_fields; $i++) {
			push @family_array, $i-9;
			# print $ped_output "$i\t";
		}
		# print $ped_output "\n";

		#individual info - use the original #CHROM header line to obtain pigID
		# print "individual_id\t";
		for (my $i = 9; $i <= $#vcf_fields; $i++) {
			push @individ_array, "$vcf_fields[$i]";
			# print $ped_output "$vcf_fields[$i]\t";
		}
		# print $ped_output "\n";

		#father_id - unknown, print 0 for all
		# print "father_id\t";
		for (my $i = 9; $i <= $#vcf_fields; $i++) {
			push @father_array, "0";
			# print $ped_output "0\t";
		}
		# print $ped_output "\n";

		#mother_id - unknown, print 0 for all
		# print "mother_id\t";
		for (my $i = 9; $i <= $#vcf_fields; $i++) {
			push @mother_array, "0";
			# print $ped_output "0\t";
		}
		# print $ped_output "\n";

		#gender - sort of known. Does it matter? Check with BL.
		# print "gender\t";
		for (my $i = 9; $i <= $#vcf_fields; $i++) {
			push @gender_array, "1";
			# print $ped_output "1\t";
		}
		# print $ped_output "\n";

		#disease status - sort of known. Does it matter? Check with BL.
		# print "status\t";
		for (my $i = 9; $i <= $#vcf_fields; $i++) {
			push @disease_array, "1"; 
			# print $ped_output "1\t";
		}
		# print $ped_output "\n";
	
	
		@mega_array = ([@family_array], [@individ_array], [@father_array], [@mother_array], [@gender_array], [@disease_array]);
	}

	else {
	
		
		my $ref_allele = $vcf_fields[3];
		my $alt_allele = $vcf_fields[4];
		
		#print info file for haploview
		print $info_output "$vcf_fields[2]\t$vcf_fields[1]\n";
		my (@genotype_array);
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

			# print $ped_output "$conversion{$individual_allele_one} $conversion{$individual_allele_two}\t";
			push @genotype_array, "$conversion{$individual_allele_one} $conversion{$individual_allele_two}";
		}
		push @mega_array, \@genotype_array;
		# print $ped_output "\n";
	}
}
		
my @transposed_array = transpose(\@mega_array);
my $dim1_length = scalar @transposed_array;
my $dim2_length = $#{$transposed_array[0]};
# print Dumper \@transposed_array;

for (my $i = 0; $i < $dim1_length; $i++) {
	for (my $j = 0; $j <= $dim2_length; $j++){
		print $ped_output "$transposed_array[$i][$j]\t";
	}
	print $ped_output "\n";
}