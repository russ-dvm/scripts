#!/usr/bin/perl

use strict;
use warnings;

#RFraser, May 27, 2016
#Take the output of the fisher's exact test script from popoolation two and grab the corresponding VCF entries from a VCF input file. Will append the significance of the FET to the VCF entry.

#Helpful if the FET file only contains significant results (could use sed for this - e.g. sed '/na/d' will delete lines containing na.)

open (my $fet, "<", $ARGV[0]);
open (my $vcf, "<", $ARGV[1]);

die ("\n\tUsage: perl <script> <fet file> <vcf file> \n\n\tThis script redirects to stdout. Redirect to a file if you want.\n\n") if @ARGV != 2;

chomp(my @fet_array = <$fet>);
chomp(my @vcf_array = <$vcf>);

#Print header
my @header = grep(/#CHROM/, @vcf_array);
print "$header[0]\tpopoolation_fet\n";




foreach my $vcf_line (@vcf_array) {

	#ignore header lines
	if ($vcf_line =~ /#/) {
	}

	else {

		my @vcf_fields = split(/\t/, $vcf_line);
		my @chrom_match = grep(/$vcf_fields[0]/, @fet_array);
		my @unique_match = grep(/$vcf_fields[1]/, @chrom_match);

		my $num_matches = @unique_match;
		if ($num_matches > 1) {

			print "Error! Multiple matches for $vcf_fields[0]\t$vcf_fields[1]\n";
		
		}
		
		elsif ($num_matches == 0) {

			print "Error! No matches for $vcf_fields[0]\t$vcf_fields[1]\n";
		}

		else {
			my $field_count = @vcf_fields;

			my @match_fields = split(/\t/, $unique_match[0]);
#			my @fet = split(/=/, $match_fields[5]);
			my $fet = $match_fields[6];
			for (my $i=0; $i<=$field_count; $i++) {
				print "$vcf_fields[$i]\t";
			}

#			print "$fet[1]";
			print "$fet";
			print "\n";
		}
}

}
