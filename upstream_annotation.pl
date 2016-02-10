#!/usr/bin/perl

use strict;
use warnings;

#this script description


open (my $meme, "<", $ARGV[0]);
open (my $output, ">", $ARGV[1]);

die "\nOops. You forgot something.\n\nProper usage: perl <script> <meme.xml file> <output>\n\n" if @ARGV != 2;


my (%genomes, %significant_motifs, @motif_id, @motif_line, @width_field, @evalue);

print "What species (horse or cow)?: ";
my $species = <STDIN>;
chomp $species;

while (my $line = <$meme>) {
	chomp $line;
	
	#Get motif name and e-value
	for (my $i = 1; $i < 16; $i++) {

		if ($line =~ m/motif id=\"motif_$i\"/) {
		
			@motif_line = split(" ", $line);
			@motif_id = split("\"", $motif_line[1]);
			@evalue = split("\"", $motif_line[8]);
			@width_field = (split"\"", $motif_line[3]);
# 			if ($evalue[1] < 0.05) {
# 				$significant_motifs{"$motif_id[1]"} = $evalue[1];
# 				print $output "$motif_id[1]\tsig\t$width_field[1]\t";
# 			}
# 			else {
# 				print $output "$motif_id[1]\tnon-sig\t$width_field[1]\t";
# 			}
		} 
	}	
	
	#Determine sequence numbers that correlate with the diff genomes
	if ($line =~ /<sequence id=\"/) {
		my @sequence_line = split("\"", $line);
		my @genome_id = split(":", $sequence_line[3]);

		#Hash - key = genome name, value equals the sequence number.		
		$genomes{"$genome_id[1]"} = "$sequence_line[1]";
	}

# 	Find positions of the bovine and equine sequences
	if ($line =~ /<contributing_site sequence_id/) {
		my @position_line = split("\"", $line);
		my $direction = $position_line[5];
		
		if ("$species" eq "horse") {
			if ("$position_line[1]" eq "$genomes{EquCab2}") {
				if ($evalue[1] < 0.05) {
# 					$significant_motifs{"$motif_id[1]"} = $evalue[1];
					print $output "$motif_id[1]\tsig\t$width_field[1]\t";
				}
				else {
					print $output "$motif_id[1]\tnon-sig\t$width_field[1]\t";
				}
				print $output "horse\t$position_line[3]\t$direction\n";
			}
		}
		
		elsif ("$species" eq "cow") {
			if ("$position_line[1]" eq "$genomes{'UMD3.1'}") {
				if ($evalue[1] < 0.05) {
# 					$significant_motifs{"$motif_id[1]"} = $evalue[1];
					print $output "$motif_id[1]\tsig\t$width_field[1]\t";
				}
				else {
					print $output "$motif_id[1]\tnon-sig\t$width_field[1]\t";
				}

				print $output "cow\t$position_line[3]\t$direction\n";
			}
		}
		else {
			print "You didn't enter a known species\n";
		}
	}
}

close $meme;
close $output;