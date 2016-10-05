#!/usr/bin/perl

use strict;
use warnings;

open (my $blocks, "<", $ARGV[0]);

my @file_array = <$blocks>;
my %conversion = (
	0 => "N",
	1 => "A",
	2 => "C", 
	3 => "G", 
	4 => "T"
	);

foreach my $line (@file_array) {

	chomp $line;
	if ($line =~ /BLOCK/) {
		print "$line\n";
	}
	else {
		my @line_fields = split(/\s/, $line);
		my @sequence = split(//, $line_fields[0]);
		foreach my $base (@sequence) {
			print "$conversion{$base}";
		}
		print " $line_fields[1]\n";
	}
}