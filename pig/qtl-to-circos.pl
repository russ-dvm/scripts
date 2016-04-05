#!/usr/bin/perl

use strict;
use warnings;

#Take eQTL results and re-format them for circos.
#requires that the raw eQTL output has been passed through the qtl-annotate.pl script to add in genomic coordinates, etc.

open (my $qtl, "<", $ARGV[0]);
open (my $link_out, ">", $ARGV[1]);
open (my $genename_out, ">", $ARGV[2]);


my @qtl_array = <$qtl>;

foreach my $qtl_line (@qtl_array) {
	chomp $qtl_line;

	next if $qtl_line =~ /snp_id/;

	my @qtl_field = split("\t", $qtl_line);

	$qtl_field[1] =~ s/chr/ss/;
	$qtl_field[4] =~ s/chr/ss/;

	my $snp_coord_adj = $qtl_field[2] + 1;

	#print LINKS - i.e. snp location to target location
	print $link_out "$qtl_field[1] $qtl_field[2] $snp_coord_adj $qtl_field[4] $qtl_field[5] $qtl_field[6]\n";

	#print GENE NAMES and their COORDINATES
	print $genename_out "$qtl_field[4] $qtl_field[5] $qtl_field[6] $qtl_field[3]\n"
}

close $qtl;
close $link_out;
close $genename_out;