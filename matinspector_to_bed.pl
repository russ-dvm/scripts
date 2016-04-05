#!/usr/bin/perl

use strict;
use warnings;


open (my $mat, "<", $ARGV[0]);
open (my $output, ">", $ARGV[1]);

die "\n\nImproper usage\n\n" if @ARGV < 2;

my $genome_dir;

while (my $mat_line = <$mat>) {
	chomp $mat_line;
	
	if ($mat_line =~ /Seq/){
	}
	else {
	
		my @mat_field = split("\t", $mat_line);
	
		#get transcription factor binding site coordinates
		my $tf_start = $mat_field[10];
		my $tf_end = $mat_field[11];
		my $tf_dir = $mat_field[13];
	
		#get genomic coordinates
		my $seq_name = $mat_field[0];
		my @seq_field = split(":", $seq_name);
			my $chrom = $seq_field[2];
			my $genome_start = $seq_field[3];
			my $genome_end = $seq_field[4];
			$genome_dir = $seq_field[5];
			my ($cis_trans);
		
	# 	create bed file
		if ($genome_dir > 0) {
			my $bed_start = $genome_start + $tf_start;
			my $bed_end = $genome_start + $tf_end;
		
			if ($tf_dir eq "+")	{
				$cis_trans = "cis";
			}
			else {
				$cis_trans = "trans";
			}
			print $output "$chrom\t$bed_start\t$bed_end\t$cis_trans\tfamily:$mat_field[4]\t$mat_field[19]\tcore_sim:$mat_field[14]\tmatrix_sim:$mat_field[15]\tinfo:$mat_field[5]\n";
		}
		
		
		elsif ($genome_dir < 0) {
			my $bed_start = $genome_end - $tf_end;
			my $bed_end = $genome_end - $tf_start;
			if ($tf_dir eq "+")	{
				$cis_trans = "cis";
			}
			else {
				$cis_trans = "trans";
			}
			print $output "$chrom\t$bed_start\t$bed_end\t$cis_trans\tfamily:$mat_field[4]\t$mat_field[19]\tcore_sim:$mat_field[14]\tmatrix_sim:$mat_field[15]\tinfo:$mat_field[5]\n";

		}
	}	
}

close $mat;
close $output;