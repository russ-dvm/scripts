#!/usr/bin/perl

use strict;
use warnings;


foreach my $file (@ARGV) {
	open (my $file_name, "<", $file);
	my @file_array = <$file_name>;
	my @gene_name_array = split(/\./, $file);
	my $printable = 0;

	foreach my $line (@file_array) {
		chomp $line;
		if ($line =~ /Alleles Captured/){
			$printable = 1;
			# print "$line\n";
		}
		elsif ($line =~ /#/){
			$printable = 0;
		}
		elsif ($line =~ /Best/) {
			$printable = 0;
		}
		elsif ($printable){
			print "$gene_name_array[0]\t$line\n";
		}
	}
}
