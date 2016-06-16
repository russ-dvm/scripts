#!/usr/bin/perl
use strict;
use warnings;

#RSF based of a script by AM, heavily modified at this point.
#rfrase03@uoguelph.ca
#2016/01/06
#This script will go through a vcf file
#and sum the reference reads and alternate reads for each of the two populations defined below.

##Modified June 13 to accept a list of normal vs diseased


my $usage = "perl script.pl vcf_file out_file\n";


open (my $vcf_file, "<", $ARGV[0]);
open (my $assignment_list, "<", $ARGV[1]);


my @vcf_array = <$vcf_file>;
my @assignments = <$assignment_list>;

#Create empty container variables
my (@header_fields);
my ($n_ref_bin, $n_alt_bin, $d_ref_bin, $d_alt_bin);


print "chrom\tposition\trsid\tref\talt\tqual\tfilter\tn.ref\tn.alt\td.ref\td.alt\n";

foreach my $line (@vcf_array) {
	chomp $line;

	##clear bins for each line
	my $n_ref_bin = 0;
	my $n_alt_bin = 0;
	my $d_ref_bin = 0;
	my $d_alt_bin = 0;


		
	if ($line =~ /\#/) {

		@header_fields = split(/\t/, $line);

	}


	else {
		my @line_fields = split("\t", $line);   	#split each entry into an array			
		for (my $i = 0; $i < scalar @line_fields; $i++) {			

			#print first 7 fields, i.e. standard VCF fields (cut out verbose info field and format)			
			
			if ($i < 7) {
				print "$line_fields[$i]\t";
			}

			elsif ($i > 8) {

				#first, find out if the genotype entry is from a diseased or normal animal
				my @dis_or_norm_array = grep { /\b$header_fields[$i]\b/ } @assignments; 
				# print "$header_fields[$i]\t";

				my @dis_or_norm = split(/\t/, $dis_or_norm_array[0]);
				if ($dis_or_norm[1] =~ /normal/) {

					# print "$line_fields[$i]\t";
					my @genotype_fields = split(/\:/, $line_fields[$i]);
					# print "$genotype_fields[1]\n";
					my @ref_count = split('\,', $genotype_fields[1]);
					$n_ref_bin = $n_ref_bin + $ref_count[0];
					$n_alt_bin = $n_alt_bin + $ref_count[1];


				}

				elsif ($dis_or_norm[1] =~ /diseased/) {
					my @genotype_fields = split(/\:/, $line_fields[$i]);
					# print "diseased\n";
					my @ref_count = split('\,', $genotype_fields[1]);
					$d_ref_bin = $d_ref_bin + $ref_count[0];
					$d_alt_bin = $d_alt_bin + $ref_count[1];

				}
			}
		}

		
		print "$n_ref_bin\t$n_alt_bin\t$d_ref_bin\t$d_alt_bin\n";

	}
}		

print "Run fet.R on the results of this script. \n";