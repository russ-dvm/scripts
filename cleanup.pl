#!/usr/bin/perl

use strict; use warnings;

#R. Fraser, Apr 13, 2015.
#Used to remove fields of irrelevant information following use of "make_diploid.pl"

open (my $vcf, "<", $ARGV[0]);						#open input vcf file
open (my $output, ">", $ARGV[1]);					#open output

while (my $line = <$vcf>) {
	chomp $line;
	
	if ($line =~ /\#/) {							#print header lines
		print $output "$line\n";
		}
	
	else {
		my @fields = split("\t", $line);			#split each snp field into an array
		my $fieldnumber = @fields;
		
		for (my $i = 0; $i < @fields; $i++) {
			if ($i < 10) {
				print $output "$fields[$i]\t";		#print the first 10 columns as they are
			}
			elsif ($fields[$i] !~ m/\d*\:\d*/g) {		#check remaining columns to see if they include the ":" character surrounded by any number of digits. Print only if the column does NOT match.
				print $output "$fields[$i]\t";
				}
			}
		print $output "\n";						#newline between entries
		}
			
}	
close $vcf;
close $output;