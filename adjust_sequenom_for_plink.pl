#!/usr/bin/perl

use strict;
use warnings;

##RFraser, 2016/06/10, rfrase03@uoguelph.ca

##Adjusts the genotypes obtained from sequenom by adding whitespace between letters, and doubling "0" if there is one. Requires a mostly prepapred .ped file (easiest to do in excel)


open (my $ped, "<", $ARGV[0]);

my @ped_array = <$ped>;

foreach my $ped_line (@ped_array) {

	chomp $ped_line;

	my @record_fields = split(/\t/, $ped_line);

	for (my $i = 0; $i < scalar @record_fields; $i++) {

		if ($i < 6) {
			print "$record_fields[$i]\t";
		}

		else {
			if ($record_fields[$i] =~ /\d+/) {
				print "0\t0\t";
			}

			else {

				my @genotype_letters = split(//, $record_fields[$i]);
				print "$genotype_letters[0]\t$genotype_letters[1]\t";
			}
		}
	}

	print "\n";

}