#!/usr/bin/perl

use strict;
use warnings;

#script description

open (my $temp, "<", $ARGV[0]);
open (my $coord_file, "<", $ARGV[1]);
open (my $output, ">", $ARGV[2]);

die "\nOops. You forgot something.\n\nProper usage: perl <script> <input file from script1> <upstream coordinate file, FASTA> <output>\n\n" if @ARGV != 3;

print $output "meme_hash=(\n";

#declare variables here
my ($meme_coord_5, $meme_coord_3, $genome_5, $genome_3, $chrom, $meme_length, $final_coord_5, $final_coord_3, $direction, $strand);

while (my $coord_file_line = <$coord_file>) {
	chomp $coord_file_line;
	
	if ($coord_file_line =~ /\>/) {
		my @header_fields = split(":", $coord_file_line);
		$chrom = $header_fields[2];
		$genome_5 = $header_fields[3];
		$genome_3 = $header_fields[4];
		$direction = $header_fields[5];
	}
		
		while (my $line = <$temp>) {
			chomp $line;
			
			my @meme_field = split("\t", $line);
			if ($meme_field[5] eq "plus"){
				$strand = "cis";
			}
			else {
				$strand = "trans";
			}

			if ($meme_field[1] eq "sig") {
				$meme_coord_5 = $meme_field[4];
				$meme_length = $meme_field[2];
				$meme_coord_3 = $meme_coord_5 + $meme_length;

				if ($direction == 1) {
					$final_coord_5 = $meme_coord_5 + $genome_5;
					$final_coord_3 = $meme_coord_3 + $genome_5 - 1;
				}
				else {
					$final_coord_3 = $genome_3 - ($meme_coord_5);
					$final_coord_5 = $genome_3 - ($meme_coord_3 - 1);
				}
				print $output "\t[\"$chrom:$final_coord_5..$final_coord_3:$direction\"\]=\"$meme_field[0]\"\n";
			}
		}
}

print $output "\t)";


close $temp;
close $coord_file;
close $output;