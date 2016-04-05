#!/usr/bin/perl

#Adapt standard VCF file to genotype input required for the matrix eQTL R package. Does both SNP genotype information and chromosomal location of the individual SNP, retaining IDs and the identical ordering.
#R. Fraser, rfrase03@uoguelph.ca

#NOTE: can NOT handle multiallelic sites. Will need to use an alternate program (e.g. GATK selectVariants or VariantsToTable) to deal with those.

use strict;
use warnings;

open (my $vcf, "<", $ARGV[0]);
open (my $output, ">", $ARGV[1]);
open (my $snpsloc, ">", $ARGV[2]);

my %genotype = (
	"0/0" => '0',
	"0/1" => '1',
	"1/1" => '2',
	"./." => 'NA',
	);

while (my $vcf_line = <$vcf>) {
	chomp $vcf_line;
	
	#ignore VCF header lines
	if ($vcf_line =~/##/){
	}
	
	#create header
	elsif ($vcf_line =~ /#CHROM/) {
	
		print $output "id\t";
		
		print $snpsloc "snp\tchr\tpos\n";
		
		my @header = split("\t", $vcf_line);
		
		for (my $i=9; $i<81; $i++) {
			print $output "$header[$i]\t";
		}
		print $output "\n";
		
	}
	
	#now get the genotypes for the SNPs
	else {
		my @line_field = split("\t", $vcf_line);
		my $id = $line_field[2];
			
		#assign ID to non-rsID SNPs	
		if ($id =~ /\./) {
			$id = "ovc" . $.;
		}
		
		#print ID (rs if available, otherwise custom OVC #)
		print $output "$id\t";
		print $snpsloc "$id\t$line_field[0]\t$line_field[1]";
		
		#figure out genotypes for each pig and print
		for (my $i=9; $i<81; $i++) {
			my @genotype_field = split(":", $line_field[$i]);	
			print $output "$genotype{$genotype_field[0]}\t";
		}

		print $output "\n";
		print $snpsloc "\n";
	}
	
}	

close $vcf;
close $output;
close $snpsloc;