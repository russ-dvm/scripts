#!/usr/bin/perl

use strict; use warnings;

#R. Fraser
#Fills the 'SNP ID' field of a VCF file with the position of the snp.

open (my $vcf, "<", $ARGV[0]);
open (my $output, ">", $ARGV[1]);

die "Usage: perl assign.snp.pl <input> <output>" unless @ARGV == 2;

while (my $line = <$vcf>) {
	chomp $line;


	if ($line =~ /\#/) {
		print $output "$line\n";					#print intro lines
		}
		
	else {
	
		
			my @fields = split("\t", $line);   				#split each entry into an array

			for (my $i = 0; $i < @fields; $i++) {			
				if ($fields[2] !~ "rs*" && $i == 2) {			#if on the third time through the loop the third column doesn't have an rsID, then
					print $output "ovc" . $. . "\t";			#change the entry (SNP ID field) to be "ovc<snp #>"
					}
				else {
					print $output "$fields[$i]\t";				#if there already is a dbSNP entry, leave it alone.
					}
				}
			print $output "\n";
			}		

}	

close $vcf;
close $output;