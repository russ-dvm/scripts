#!/usr/bin/perl

use strict;
use warnings;

#R.Fraser, rfrase03@uoguelph.ca, 2016/02/05. 
#Takes input from Promo3, isolates factors with a predicted RE (query) value of less than a prompted value, and then outputs in either list or bed format.

open (my $tf_file, "<", $ARGV[0]);
open (my $output, ">", $ARGV[1]);

die "\nImproper usage.\n\nTry: perl <script> <input_file_from_promo> <output>\n\n" if @ARGV != 2;


##failed attempt to use getopt::long

	# my $format = "bed";
	# my $cutoff = 100;
	# my ($tf_file, $output);
	# 
	# GetOptions(
	# 	"format=s" => \$format,
	# 	"cutoff=i" => \$cutoff,
	# 	"input=s" => \$tf_file,
	# 	"output=s" => \$output)
	# or die ("\nUsage incorrect.\n
	# perl make_tf_sites_bed.pl --input <input> --output <output>\n
	# Required:
	# \t --input\tInput file from Promo3
	# \t --output\tOutput file
	# 
	# Options:
	# \t--format\tSpecify if you want in list or bed format
	# \t\t\t(Default: bed)
	# \t--cutoff\tSpecify the cutoff value of RE
	# \t\t\t(Default: 100)\n\n");
	# 



	

print "Format (list or bed)?: ";
my $format = <STDIN>;
chomp $format;

print "Cutoff value?: ";
my $cutoff = <STDIN>;
chomp $cutoff;

if (($format eq "list") or ($format eq "bed")) {
	while (my $line = <$tf_file>) {
	
	
		chomp $line;	
		my @field = split(";", $line);
	
		if ($line =~ /chromo/) {
	
			my @genome_info = split(":", $field[0]);
				my $chrom = $genome_info[2];
				my $start_pos = $genome_info[3];
	
			my $trans_factor = $field[1];	
	
			my $trans_factor_start = $field[2];
				$trans_factor_start =~ s/^\s+//;
			my $trans_factor_end = $field[3];
				$trans_factor_end =~ s/^\s+//;
			my $trans_factor_seq = $field[5];

			my $re_query = $field[7];
			

			if ($re_query < 5){

				my $genome_start = "$start_pos" + "$trans_factor_start";
				my $genome_end = "$start_pos" + "$trans_factor_end";

				#list format		
				if ($format eq "list") {
					print $output "$chrom:$genome_start-$genome_end\t$trans_factor |$trans_factor_seq |$re_query\n";
				}

				#bed format
				elsif ($format eq "bed") {
					print $output "chr$chrom\t$genome_start\t$genome_end\t$trans_factor |$trans_factor_seq |$re_query\n";
				}
			}
		}
	}
}

else {
	die "\nWrong value entered. Bye.\n\n";
	}

close $tf_file;
close $output;