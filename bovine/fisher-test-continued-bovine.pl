#!/usr/bin/perl
use strict;
use warnings;

#RSF based of a script by AM.
#rfrase03@uoguelph.ca
#2016/01/06
#This script will go through a vcf file
#and sum the reference reads and alternate reads for each of the two populations defined below.


my $usage = "perl script.pl vcf_file out_file\n";

my $vcf_file = shift || die $usage;
my $out_file = shift || die $usage;


open IN, "$vcf_file" || die "could not open $vcf_file for reading\n";
open OUT, ">$out_file" || die $usage;
print OUT "chrom\tposition\trsid\tref\talt\tqual\tfilter\tn.ref\tn.alt\td.ref\td.alt\n";

while (<IN>) {
	s/\cM//;
	chomp;

	my $line = $_;
	my $n_ref_bin = 0;
	my $n_alt_bin = 0;
	my $d_ref_bin = 0;
	my $d_alt_bin = 0;

		
	if ($line =~ /\#/) {
		}
		
	else {
		my @fields = split("\t", $line);   	#split each entry into an array			
		for (my $i = 0; $i < 7; $i++) {						
			print OUT "$fields[$i]\t";		#print first 7 fields, i.e. standard VCF fields (cut out verbose info field and format)
		}
		my $count = $fields[9];		#variable count is assigned the value of the string of read counts
		my @depth = split('\"', $count);				
		for (my $i=0; $i<21; $i++) {
			if ($i==0){		#NORMAL, group1
				my @ref_count = split('\,', $depth[$i]);
				$n_ref_bin = $n_ref_bin + $ref_count[0];
				$n_alt_bin = $n_alt_bin + $ref_count[1];
				}
			elsif ($i>0 && $i<10) {		#DISEASED, groups11-18
				my @ref_count = split('\,', $depth[$i]);
# 				print "$ref_count[0]\t";
# 				print "$ref_count[1]\n";
				$d_ref_bin = $d_ref_bin + $ref_count[0];
				$d_alt_bin = $d_alt_bin + $ref_count[1];
				}
			elsif ($i>9 && $i<12) {			#NORMAL, group19, 2
				my @ref_count = split('\,', $depth[$i]);
# 				print "$ref_count[0]\t";
# 				print "$ref_count[1]\n";
				$n_ref_bin = $n_ref_bin + $ref_count[0];
				$n_alt_bin = $n_alt_bin + $ref_count[1];
				}
			elsif ($i>12 && $i <18) {			#NORMAL, groups3-6, 24
				my @ref_count = split('\,', $depth[$i]);
# 				print "$ref_count[0]\t";
# 				print "$ref_count[1]\n";
				$n_ref_bin = $n_ref_bin + $ref_count[0];
				$n_alt_bin = $n_alt_bin + $ref_count[1];
				}
			else {					#DISEASED, remainder
				my @ref_count = split('\,', $depth[$i]);
# 				print "$ref_count[0]\t";
# 				print "$ref_count[1]\n";
				$d_ref_bin = $d_ref_bin + $ref_count[0];
				$d_alt_bin = $d_alt_bin + $ref_count[1];
				}
			}
		
		print OUT "$n_ref_bin\t$n_alt_bin\t$d_ref_bin\t$d_alt_bin\n";

	}
}		

print "Run fet.R on the results of this script. \n";

close IN;
close OUT; 