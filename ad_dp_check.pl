#!/usr/bin/perl
use strict;
use warnings;

#RSF based of a script by AM.
#rfrase03@uoguelph.ca
#2016/01/26
#This script checks to see if the sum of the allele depths for a specified group manages its reported DP, as per the issues described in http://gatkforums.broadinstitute.org/gatk/discussion/6005/allele-depth-ad-is-lower-than-expected

my $usage = "perl script.pl vcf_file out_file\n";

my $vcf_file = shift || die $usage;
my $out_file = shift || die $usage;

open IN, "$vcf_file" || die "could not open $vcf_file for reading\n";
open OUT, ">$out_file" || die $usage;

print OUT "chrom\tpos\tnum.ref\tnum.alt\tgatk.total\tcalc.total\tequal\n";


while (<IN>) {
	s/\cM//;
	chomp;

	my $line = $_;
	my $total = 0;
	if ($line =~ m/^#/) {
	} 
		
	else {
		my @fields = split("\t", $line);   	#split each entry into an array			
		for (my $i = 0; $i < 2; $i++) {						
			print OUT "$fields[$i]\t";		#print first 9 fields, i.e. standard VCF fields
			}
		for (my $i=30) {	#sloppy carry-over from another script... just change i so that it only equals one number
			my $group_data = $fields[$i];		
			my @split_data = split('\:', $group_data);
			my @reads = split('\,', $split_data[1]);
			$total = $reads[0] + $reads[1];
			print OUT "$reads[0]\t$reads[1]\t";
			print OUT "$split_data[2]\t";
			print OUT "$total\t";
			if ($total = $split_data[2]){
				print OUT "TRUE\n";
				}
			else {
				print OUT "FALSE\n";
				}
			}
	}
}		

close IN;
close OUT; 	
	
print "Next, run fisher-test-continued.pl on the output of this program\n";