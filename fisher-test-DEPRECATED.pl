#!/usr/bin/perl
use strict;
use warnings;

#RSF based of a script by AM.
#rfrase03@uoguelph.ca
#2016/01/06
#This script will go through a VCF file and remove the genotype information, concatenating
# reads for each group and separating them with a ".

my $usage = "perl script.pl vcf_file out_file\n";

my $vcf_file = shift || die $usage;
my $out_file = shift || die $usage;

open IN, "$vcf_file" || die "could not open $vcf_file for reading\n";
open OUT, ">$out_file" || die $usage;
while (<IN>) {
	s/\cM//;
	chomp;

	my $line = $_;

	if ($line =~ m/^#/) {
		print OUT "$line\n";
	} 
		
	else {
		my @fields = split("\t", $line);   	#split each entry into an array			
		for (my $i = 0; $i < 9; $i++) {						
			print OUT "$fields[$i]\t";		#print first 9 fields, i.e. standard VCF fields
			}
		for (my $i=9; $i<30; $i++) {
			my $group_data = $fields[$i];		
			my @split_data = split('\:', $group_data);
			my @reads = split('\,', $split_data[1]);
			print OUT "$reads[0],$reads[1]\"";
			}
		print OUT "\n";
	}
}		

close IN;
close OUT; 	
	
print "Next, run fisher-test-continued.pl on the output of this program\n";