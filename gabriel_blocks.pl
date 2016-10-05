#!/usr/bin/perl

use strict;
use warnings;

open (my $gabriel, "<", $ARGV[0]);

my @gabriel_array = <$gabriel>;

my %translate_hash = (
	"1" => "A",
	"2" => "C",
	"3" => "G",
	"4" => "T");

foreach my $line (@gabriel_array){
	chomp $line;
	if ($line =~ /BLOCK/){}
	else {
		my @full_sequence = split(/ /, $line);
		my @sequence = split(//, $full_sequence[0]);
		foreach my $letter (@sequence){
			print "$translate_hash{$letter}"
		}
		print "\t$full_sequence[1]\n";
	}
}