#!/bin/usr/perl

use strict;
use warnings;

# Annotates the matrix eQTL output with SNP locations and Gene locations.

open (my $qtl, "<", $ARGV[0]);
open (my $snpsloc, "<", $ARGV[1]);
open (my $geneloc, "<", $ARGV[2]);
open (my $output, ">", $ARGV[3]);

die "\n Usage: perl <script> <qtl_file> <snps_loc_file> <gene_loc_file> <output> \n\n" if @ARGV != 4;


my @snpsloc_lines = <$snpsloc>;
my @geneloc_lines = <$geneloc>;



while (my $qtl_line = <$qtl>) {
	
	chomp $qtl_line;
	
	#fix header line here
	if ($qtl_line =~ /snps/) {
		print $output "snp_id\tsnp_chrom\tsnp_position\tgene\tgene_chrom\tgene_start\tgene_end\tstatistic\tpvalue\tFDR\tbeta\n";
	}
	else {
		my @qtl_fields = split("\t", $qtl_line);
		
		#get snp id from the QTL file	
		my $snp_id = $qtl_fields[0];
		
		#find the location info from the snpsloc file
		my @snp_location = grep(/$snp_id/, @snpsloc_lines);
		chomp $snp_location[0];
		
		#print snp location info to output	
		print $output "$snp_location[0]\t";
		
		#get gene id
		my $gene_id = $qtl_fields[1];
		
		#find gene location info from geneloc file
		my @gene_location = grep(/$gene_id/, @geneloc_lines);

		##print gene location info to output
		#test to see if the gene was discovered
		if (@gene_location) {
		chomp $gene_location[0];
		print $output "$gene_location[0]\t";
		}
		else {
		print $output "$gene_id\tNA\tNA\tNA\t";
		}
		
		#print remainder of the QTL info
		print $output "$qtl_fields[2]\t$qtl_fields[3]\t$qtl_fields[4]\t$qtl_fields[5]\n";	
		
	}
}

