#!/usr/bin/perl

use strict;
use warnings;
use Statistics::Basic qw(:all); # this package calculates population standard deviation rather than sample deviation. To go back to using it change the calcs to be stddev(array).


#RFraser, rfrase03@uoguelph.ca
#April 03, 2016
#This script should do the following:
#	1) Access the results of Matrix eQTL to identify significant SNP-expression associations. 
# 	***NOTE: USE THE OUTPUT OF QTL-ANNOTATE.pl for this script***
#	2) Identify the SNP_id and gene_id of the significant associations
#	3) Take the SNP_id and harvest the genotypes from all animals for that id from the genotype data derived from the Illumina sequencing run (a file noted as SNP-(date).txt).
#	4) Take the gene_id and harvest the expression data for all the animals from the microarray data file prepared for Matrix eQTL (shoudl be found a "gene_expression_(data).txt").
#	5) Form a table consisting of 3 columns: animal_id	genotype	expression
#	6) Summarize the above data based on genotype: provide another table with total number of animals genotyped, how many of each genotype (homozygous ref/alt or hetero), avg expression, stddev of expression, ref and alt allele, position of the gene and of the snp.
#	7) For the sake of simplicity, generate a seperate output file for each SNP_id (name is based on the SNP_id). If all the results are desired in one file can write a simple bash script to accomplish this later. These individual files will be useful when importing into R and computing boxplots etc of genotype vs expression.


open (my $genotypes, "<", $ARGV[0]);
open (my $expression, "<", $ARGV[1]);
open (my $qtls, "<", $ARGV[2]);
open (my $vcf, "<", $ARGV[3]);

die "\nUsage: perl <script> genotype.file(SNPs) expressionFile(expression_data_date) qtlfile(output from eQTL) vcf_file(vcf of variants\n\n" if @ARGV != 4;
#add a die statement and maybe a usage statement... god knows i'm goign to forget what's required

my @qtl_array = <$qtls>;
my @genotype_array = <$genotypes>;
my @expression_array = <$expression>;
my @vcf_array = <$vcf>;

my @genotype_header = split("\t", $genotype_array[0]);
my @expression_header = split("\t", $expression_array[0]);

my $shitty_solution = 0;


##need to fix the counter to reflect $#genotype_array
#start at one to skip the header line
for (my $i=1; $i <= $#qtl_array; $i++) {

	#generate a counter to uniquify files later on. The problem is that a certain number of SNPs affect more than one locus; they need to be uniqufied in order to prevent overwriting earlier iterations. On the plus side, this number also represents the rank of the SNP in terms of significance, so not totally awful I guess...
	$shitty_solution = $shitty_solution + 1; 

	#set counters and arrays to zero for each iteration through the files
	my $reference_total = 0;
	my $reference_counter = 0;
	my @reference_array = ();

	my $hetero_total = 0;
	my $hetero_counter = 0;
	my @hetero_array = ();

	my $alternate_total = 0;
	my $alternate_counter = 0;
	my @alternate_array = ();

	my $na_count = 0;

	#Identify SNP_ids and their matching genes from the eQTL output file
	my @qtl_fields = split("\t", $qtl_array[$i]);
	my $snp_of_interest = $qtl_fields[0];
	my $gene_of_interest = $qtl_fields[3];

	#Find the genotype for each SNP
	my @matched_genotype_line = grep {/$snp_of_interest/} @genotype_array;
	#Find the expression data for each gene associated with the above snp
	my @matched_expression_line = grep {/$gene_of_interest/} @expression_array;


	#Need a method of handling SNP that have matched more than one gene.
	my $num_of_snp_matches = grep {/$snp_of_interest/} @qtl_array;
	#num_of_snp_matches reflects the number of time that the snp appears in the qtl results.


	my @genotype_data = split("\t", $matched_genotype_line[0]);
	my @expression_data = split("\t", $matched_expression_line[0]);


	#transpose the genotypes and expression data into columns

	# #create output file name based on snpID
	my $outfile = "$genotype_data[0].$shitty_solution.txt";
	open (my $output, ">", "files/$outfile");

	#Header line for the genotype & expression data table
	print $output "animal\tgenotype\texpression.$expression_data[0]\tgenotype_letters\n";
	
	#get the reference and alternate allele from the original VCF file
	my @vcf_match = grep {/$qtl_fields[2]/} @vcf_array;
	my @vcf_fields = split("\t", $vcf_match[0]);
	my $ref_allele = $vcf_fields[3];
	my $alt_allele = $vcf_fields[4];
	my ($genotype_letters);

	#Print the individual pig genotype and expression into rows
	for (my $i=1; $i <= $#genotype_data; $i++) {
		chomp $genotype_header[$i];
		chomp $genotype_data[$i];
		chomp $expression_data[$i];
		
		#set up bins and counters for averages and stdevs, and totals animals genotyped
		if ($genotype_data[$i] eq "NA") {
			$na_count = $na_count + 1;
			$genotype_letters = "NA";
		}
		
		elsif ($genotype_data[$i] == 0) {
			$reference_total = $reference_total + $expression_data[$i];
			$reference_counter = $reference_counter + 1;
			push @reference_array, "$expression_data[$i]";
			$genotype_letters = "$ref_allele$ref_allele";
		}
		elsif ($genotype_data[$i] == 1) {
			$hetero_total = $hetero_total + $expression_data[$i];
			$hetero_counter = $hetero_counter + 1;
			push @hetero_array, "$expression_data[$i]";
			$genotype_letters = "$ref_allele$alt_allele";
		}
		elsif ($genotype_data[$i] == 2) {
			$alternate_total = $alternate_total + $expression_data[$i];
			$alternate_counter = $alternate_counter + 1;
			push @alternate_array, "$expression_data[$i]";
			$genotype_letters = "$alt_allele$alt_allele";
		}

		#print data (by row) - genotype-header = pig names.
		print $output "$genotype_header[$i]\t$genotype_data[$i]\t$expression_data[$i]\t$genotype_letters\n";
	}

	#calculate averages & stdevs for diff genotypes
	my $rr_avg = mean(@reference_array);
	my $rr_stdev = stddev(@reference_array);

	my $ra_avg = mean(@hetero_array);
	my $ra_stdev = stddev(@hetero_array);

	my $aa_avg = mean(@alternate_array);
	my $aa_stdev = stddev(@alternate_array);

	my $total_genotyped = $reference_counter + $hetero_counter + $alternate_counter;

	# #get the reference and alternate allele from the original VCF file
	# my @vcf_match = grep {/$qtl_fields[2]/} @vcf_array;
	# my @vcf_fields = split("\t", $vcf_match[0]);
	# my $ref_allele = $vcf_fields[3];
	# my $alt_allele = $vcf_fields[4];


	#set up second summary table
	print $output "\n";

	#print header line of summary table
	print $output "#gene\tgene_chrom\ttranscript_start\ttranscript_end\tsnp_id\tsnp_chrom\tsnp_pos\tref_allele\talt_allele\tn_rr\tavg_rr\tstdev_rr\tn_ra\tavg_ra\tstdev_ra\tn_aa\tavg_aa\tstdev_aa\ttotal_num_genotyped\ttotal_na\tp-value\tFDR\tslope\n";
	
	#print data
	print $output "$expression_data[0]\t$qtl_fields[4]\t$qtl_fields[5]\t$qtl_fields[6]\t$genotype_data[0]\t$qtl_fields[1]\t$qtl_fields[2]\t$ref_allele\t$alt_allele\t$reference_counter\t$rr_avg\t$rr_stdev\t$hetero_counter\t$ra_avg\t$ra_stdev\t$alternate_counter\t$aa_avg\t$aa_stdev\t$total_genotyped\t$na_count\t$qtl_fields[8]\t$qtl_fields[9]\t$qtl_fields[10]";
	
}	
