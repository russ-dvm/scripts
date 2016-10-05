#!/usr/bin/perl

use strict;
use warnings;
use List::MoreUtils "uniq";


open(my $vcf, "<", $ARGV[0]);

my @vcf_array = <$vcf>;
chomp @vcf_array;
my (@gene_names);

my @collectins = qw(MBL1 MBL2 FCN1 FCN2 FCN3 SFTPA1 SFTPD COLEC10 COLEC11 COLEC12 MASP1 MASP2);
my $file_length = @vcf_array;

#skip header line, $i = 1
for (my $i = 1; $i < $file_length; $i++){

	my @vcf_fields = split("\t", $vcf_array[$i]);
	my @gene_probes = split(/\./, $vcf_fields[0]);

	push @gene_names, $gene_probes[0];
}

my @merged = (@gene_names,@collectins);

my @unique_merged = uniq @merged;


foreach my $gene (@merged) {
	my @search_result = grep(/$gene/, @vcf_array);



	if (@search_result){
		my $filehandle = "$gene.list";
		open(my $output, ">", $filehandle);

		foreach my $line (@search_result) {
			my @line_field = split(/\t/, $line);
			print $output "$line_field[4]\n";
	
		}
	}
}