#!/usr/bin/perl

use strict;
use warnings;

open (my $ref, "<", $ARGV[0]);
open (my $fixer, "<", $ARGV[1]);

my @ref_array = <$ref>;
my @fixer = <$fixer>;
my $gene;
foreach my $fixer_line (@fixer) {

	chomp $fixer_line;
	#ignore header line
	if ($fixer_line =~ /chrom/) {}
	
	#parse the fixer file
	else {
		my @fixer_line_array = split(/\t/, $fixer_line);
		my $fixer_Chrom = $fixer_line_array[0];
		my $fixer_Start = $fixer_line_array[1];
		my $fixer_End = $fixer_line_array[2];


		#find entries that match the chrom in the file
		my @chrom_matches = grep(/$fixer_Chrom/, @ref_array);

		foreach my $chrom_match_line (@chrom_matches){
			chomp $chrom_match_line;

			my @chrom_match_line_fields = split(/\t/, $chrom_match_line);
			#determine strand so that coordinates can be properly searched
			if ($chrom_match_line_fields[3] == 1) {

				my $tss_Start = $chrom_match_line_fields[8];
				my $tss_End = $chrom_match_line_fields[9];

				if (($fixer_End > $tss_Start and $fixer_End < $tss_End) or ($fixer_Start > $tss_Start and $fixer_Start < $tss_End)){
					# print "$chrom_match_line\n"
					$gene = $chrom_match_line_fields[5];
					# print "$gene\n";
				}

			}

			elsif ($chrom_match_line_fields[3] < 0) {
				my $tss_Start = $chrom_match_line_fields[8];
				my $tss_End = $chrom_match_line_fields[9];

				if (($fixer_End > $tss_Start and $fixer_End < $tss_End) or ($fixer_Start > $tss_Start and $fixer_Start < $tss_End)){
					$gene = $chrom_match_line_fields[5];
					# print "$gene\n"

				}
			}
		}
		print "$fixer_line$gene\n";
	}
}
