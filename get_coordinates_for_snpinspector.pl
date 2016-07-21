#!/usr/bin/perl

use strict;
use warnings;

open (my $vcf, "<", $ARGV[0]);
open (my $coord, "<", $ARGV[1]);

my %rev_comp = (
	'A' => 'T',
	'T' => 'A',
	'C' => 'G',
	'G' => 'C',
	);

my ($genomatix_start, $genomatix_end, $genomatix_dir, $genomatix_chrom, $allelePos, $allelePosAdjust, $gene_name, @sequence);
my ($vcf_chrom, $var_locus, $ref_allele, $alt_allele);

while (my $geno_line = <$coord>){
	chomp $geno_line;
	
	if ($geno_line =~ /\>/) {
		#get info
		my @geno_field = split(/\|/, $geno_line);
		my @pre_geno_chrom = split("=", $geno_field[6]);
			$genomatix_chrom = "$pre_geno_chrom[0]$pre_geno_chrom[1]";
		my @pre_geno_start = split("=", $geno_field[9]);
			$genomatix_start=$pre_geno_start[1];
		my @pre_geno_end = split("=", $geno_field[10]);
			$genomatix_end=$pre_geno_end[1];
		my @pre_gene_name = split("=", $geno_field[1]);
			$gene_name = $pre_gene_name[1];
		my @pre_dir = split("=", $geno_field[8]);
			$genomatix_dir = $pre_dir[1];
		$allelePosAdjust = $genomatix_start - 1;

	}
	#store the remainder of the sequence for use later.
	else {
		push @sequence, $geno_line;
	}	
}

while (my $line = <$vcf>){
	chomp $line;

	if ($line =~ /\#/) {
	}
	else {
		#get the info
		my @vcf_fields = split("\t", $line);
		$vcf_chrom = $vcf_fields[0];
		$var_locus = $vcf_fields[1];
		$ref_allele = $vcf_fields[3];
		$alt_allele = $vcf_fields[4];
		
		#match the variants to the genomatix file
		if (($vcf_chrom eq $genomatix_chrom) && ($var_locus < $genomatix_end) && ($var_locus > $genomatix_start)) {
		
			#determine allele position using the start of the genomatix file as base 1
			$allelePos = $var_locus - $allelePosAdjust;
		
			#account for reverse complement
			if ($genomatix_dir eq "(-)"){
				$allelePos = $genomatix_end - $var_locus + 1;
				my $a = length($ref_allele);
				my $b = length($alt_allele);

				#generate reverse complement for insertions/deletions
				#deletions
				if ($a > 1) {
					my @bases = split(//, $ref_allele);
					my @ref_base_container;
					print "$ref_allele\t";
					foreach my $base (@bases) {
						unshift @ref_base_container, $rev_comp{$base};
						$ref_allele = "@ref_base_container";
						$ref_allele =~ s/\s//g;
						#adjust allelePos beceause insertions push things back
						$allelePos = ($genomatix_end - $var_locus + 1) - ($a - 1);
						
					}
# 					print "@ref_base_container\n";
					print "$ref_allele\n";
				}
				elsif ($a == 1) {
					$ref_allele = $rev_comp{$ref_allele};
				}

				#insertions, alternate alleles
				if ($b > 1) {
					my @bases = split(//, $alt_allele);
					my @alt_base_container;
					foreach my $base (@bases) {
						unshift @alt_base_container, $rev_comp{$base};
						$alt_allele = "@alt_base_container";
						$alt_allele =~ s/\s//g;
					}
				}
				
				#regular alt alleles
				elsif ($b == 1) {
 					$alt_allele = $rev_comp{$alt_allele};
				}
			
		}	
			
		#print a buncha files
		my $title = "$gene_name\_$vcf_chrom\_$allelePos.fa";
		open (my $out, '>', $title);
		print $out ">$gene_name|allelePos=$allelePos|alleles=\"$ref_allele/$alt_allele\"|$vcf_chrom|ref_pos=$var_locus\n";
		print $out @sequence;
		print $out "\n";
		}
	
	}

}

close $coord;
close $vcf;
