#!/usr/bin/perl

use strict;
use warnings;

open (my $sequenom, "<", $ARGV[0]);
open (my $list, "<", $ARGV[1]);
open (my $output, ">", $ARGV[2]);


my %bases = (
		A => "T",
		T => "A",
		C => "G",
		G => "C",
		);

my @sequenom_lines = <$sequenom>;
my @list_lines = <$list>;

print $output "$sequenom_lines[0]";

foreach my $line (@list_lines) {
	chomp $line;
	
	
	foreach my $i (1 .. $#sequenom_lines) {
		chomp $sequenom_lines[$i];

		my @sequenom_fields = split("\t", $sequenom_lines[$i]);

		if ($sequenom_fields[0] eq $line) {
			print $output "$sequenom_fields[0]\t";
			
			foreach my $i (1 .. $#sequenom_fields) {
				my @sequenom_geno = split('', $sequenom_fields[$i]);
				$sequenom_fields[$i]= "$bases{$sequenom_geno[0]}$bases{$sequenom_geno[1]}";

				print $output "$sequenom_fields[$i]\t";
			}
		}
	}
	print $output "\n";
}

close $sequenom;
close $list;
close $output;