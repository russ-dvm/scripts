#!/usr/bin/perl

use strict;
use warnings;


open (my $vcf_file, "<", $ARGV[0]);
open (my $bed_file, "<", $ARGV[1]);
open (my $output, ">", $ARGV[2]);
# 
die "\nWrong usage.\n\nperl identify_promoter_snps <vcf-input> <bed-input> <output>\n\n" if @ARGV < 3;

my ($snp_position, $start_pos, $end_pos, @snps);

while (my $vcf_line = <$vcf_file>) {

	chomp $vcf_line;
	
	if ($vcf_line =~ /\#/) {
	}
	else {
		my @vcf_fields = split("\t", $vcf_line);
		push @snps, $vcf_fields[1];
	}
}


while (my $bed_line = <$bed_file>) {
	chomp $bed_line;

	my @bed_fields = split("\t", $bed_line);
	$start_pos = $bed_fields[1];
	$end_pos = $bed_fields[2];
	
	foreach my $snp (@snps) {
		if ( $snp >= $start_pos && $snp <= $end_pos) {
			print $output "$bed_line | snp=$snp\n";
		}
	}
}

close $vcf_file;
close $bed_file;
close $output;
	
		