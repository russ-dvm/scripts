use strict;
use warnings;

#RSF, rfrase03@uoguelph.ca
#Issue: after running matrix eQTL, its apparent that if there is only ONE snp of a certain genotype, it can significantly affect the results.
#This file attempts to adjusts the genotypes in teh SNP file used by matrix eQTL such that
# 1) if there are < 3 animals in any genotype, they get set to NA
# 2) is done in a semi-intelligent manner, such that the remainder of the genotypes are left alone.

open (my $snps_with_genotypes, "<", $ARGV[0]);
open (my $snp_list_less_than_three, "<", $ARGV[1]);
open (my $output, ">", $ARGV[2]);

my @snp_with_genotype_array = <$snps_with_genotypes>;
my @snp_list_less_than_three_array = <$snp_list_less_than_three>;

foreach my $snp_output_line (@snp_with_genotype_array) {

	chomp $snp_output_line;

	#print header
	if ($snp_output_line =~ /id/){
		print "$snp_output_line\n";
	}

	else {

		###match the RSID from the SNP list to the low-genotype list\
		#first isolate the rsID from the genotyped snps
		my @snp_output_fields = split("\t", $snp_output_line);
		my $output_snp_id = $snp_output_fields[0];

		#nifty trick using the \b boundary reg-ex - insures that the exact variable is met rather than a partila match (e.g. ovc277 will only match ovc277 and not ovc27771)
		if (grep {/\b$output_snp_id\b/} @snp_list_less_than_three_array) {
			
			#these SNPs have at least one genotype that has < 3 entries

			#need to determine which genotype - dominant, hetero, recessive (these are labelled in column two of the input file)

			#store the result

			my @grep_result = grep {/\b$output_snp_id\b/} @snp_list_less_than_three_array;
			chomp @grep_result;
			my @grep_result_fields = split("\t", $grep_result[0]);

			if ($grep_result_fields[1] eq "dominant"){
				foreach my $line_field (@snp_output_fields) {
					if ($line_field eq 0) {
						print "NA\t";
					}
					else {
						print "$line_field\t"
					}
				}
				print "\n";
			}

			elsif ($grep_result_fields[1] eq "heterozygous") {
				foreach my $line_field (@snp_output_fields) {
					if ($line_field eq 1) {
						print "NA\t";
					}
					else {
						print "$line_field\t"
					}
				}
				print "\n";
			}

			elsif ($grep_result_fields[1] eq "recessive"){
				foreach my $line_field (@snp_output_fields) {
					if ($line_field eq 2) {
						print "NA\t";
					}
					else {
						print "$line_field\t"
					}
				}
				print "\n";
			}

		}

		else {

			print "$snp_output_line\n";
		}

		#requires isolating snp_ids from both input arrays

		#if doesn't match, print line

		#if does match, make modifications then print the line;

	}
}