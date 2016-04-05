#!/usr/bin/perl

use strict;
use warnings;

open (my $sequenom, "<", $ARGV[0]);
open (my $illumina, "<", $ARGV[1]);
open (my $out, ">", $ARGV[2]);

my @sequenom_lines = <$sequenom>;
my @illumina_lines = <$illumina>;

print $out "$sequenom_lines[0]";

foreach my $i (1 .. $#sequenom_lines) {
	chomp $sequenom_lines[$i];
	my @sequenom_fields = split ("\t", $sequenom_lines[$i]);
	
	foreach my $i (1 .. $#illumina_lines) {
		
		chomp $illumina_lines[$i];
		
		my @illumina_fields = split("\t", $illumina_lines[$i]);
		
		if ($illumina_fields[0] eq $sequenom_fields[0]){
			
			print $out "$illumina_fields[0]\t";
						
			foreach my $i (1 .. $#sequenom_fields) {
	
# 				deal with no-call (..) genotypes				
				if ($illumina_fields[$i] eq ".." or $sequenom_fields[$i] eq "..") {
					print $out "NA\t";
				}
			
# 				genotypes that match exactly (homozygotes, some heterozygotes)
				elsif ($sequenom_fields[$i] eq $illumina_fields[$i]) {
					print $out "1\t";
				}
			
# 				match flipped heterozygotes (e.g. A/T = T/A)
				elsif ($sequenom_fields[$i] eq reverse $illumina_fields[$i]){
					print $out "1\t";
				}
			

# 				non-matching												
				else {
					print $out "0\t";
				}
			
# 				print $out "$sequenom_fields[$i] $illumina_fields[$i]\t";
	
			}
		print $out "\n";
		}
	}	
}