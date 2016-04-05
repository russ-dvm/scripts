#!/usr/bin/perl

use strict;
use warnings;

#Quick script to add a column including the genotype letters to the files generated from "brandons_table.pl"


open (my $file, "<", $ARGV[0]);

chomp (my @file_lines = <$file>);

my ($ref_allele, $alt_allele);

for (my $i=75) {
	my @results_array = split("\t", $file_lines[$i]);
	$ref_allele = $results_array[7];
	$alt_allele = $results_array[8];
	print "$ref_allele $alt_allele\n";
}

for (my $i=0; $i=75; $i++) {
	if ($i == 0) {
		print "$file_lines[$i]\tgeno\n";
		}
}
