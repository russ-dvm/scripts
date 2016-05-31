#!/usr/bin/perl

#RSFraser, 2016/05/19
#Harvest transcript IDs (ensembl) from a SnpEFF annotated VCF File for use in batch query with Polyphen2
#This script is almost identical to transcripts_from_vcf.pl, except it includes additional info (ref/alt allele, chr, rsID, etc), and is in tab delimited format for import into excel.

open (my $input, "<", $ARGV[0]);

my @input_array = <$input>;

my %amino_acid = (
	Ala => "A",
	Arg => "R",
	Asn => "N",
	Asp => "D",
	Cys => "C",
	Glu => "E",
	Gln => "Q",
	Gly => "G",
	His => "H",
	Hyp => "O",
	Ile => "I",
	Leu => "L",
	Lys => "K",
	Met => "M",
	Phe => "F",
	Pro => "P",
	Glp => "U",
	Ser => "S",
	Thr => "T",
	Trp => "W",
	Tyr => "Y",
	Val => "V"
);

die ("\nPlease specify input file\n\n") if @ARGV != 1;

print "transcript_id\taa_position\tref_aa\talt_aa\tchrom\tsnp_pos\trsID\tref_allele\talt_allele\n";

foreach $input_line (@input_array) {
	chomp $input_line;
	my @input_field = split("\t", $input_line);
	my $annotate_field = $input_field[7];
	my @snp_eff_field = split("ANN", $annotate_field);

	my @snp_effect = split(/\|\|/, $snp_eff_field[1]);

	foreach my $effect (@snp_effect) {
		if ($effect =~ /missense_variant/) {
			my @output_info = split(/\|/, $effect);

			#clean up transcript ID
			my @transcript_id = split(/\./, $output_info[6]);
			my $transcript_id_final = $transcript_id[0];

			#reorganize the amino acids and change 3 letter code -> 1 letter
			my @aa_array = split(/\./, $output_info[10]);
			my @aa_blob = split(/(\d+)/, $aa_array[1]);

			my $ref_aa = $amino_acid{$aa_blob[0]};
			my $alt_aa = $amino_acid{$aa_blob[2]};
			my $position = $aa_blob[1];

			print "$transcript_id_final\t$position\t$ref_aa\t$alt_aa\t$input_field[0]\t$input_field[1]\t$input_field[2]\t$input_field[3]\t$input_field[4]\n";
		
		}
	}
}

