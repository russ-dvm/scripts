#!/usr/bin/perl

use strict;
use warnings;
use List::MoreUtils qw(firstidx);

#Idenitfy which AHL pigs are in which groups
#Determine the genotype of the 6-ploid organisms.
#Count 'em up just like pooling_vs_sequenom.pl.

open (my $ahl_list, "<", $ARGV[0]);
open (my $sequenom, "<", $ARGV[1]);

my @ahl_array = <$ahl_list>;
my @sequenom_array = <$sequenom>;
my @sequenom_header = split(/\t/, $sequenom_array[0]);
my $index;
my (@ahl_entry_field, @sequenom_genotypes, @genos);


foreach my $sequenom_variant (@sequenom_array){
	if ($sequenom_variant =~ /rsid/){}

	else {
		chomp $sequenom_variant;
		@sequenom_genotypes = split(/\t/, $sequenom_variant);

		foreach my $ahl_entry (@ahl_array[1..$#ahl_array]) {
			chomp $ahl_entry;

			if ($ahl_entry =~ /ID/){print "hello\t";
			}

			else {
				@ahl_entry_field = split(/\t/, $ahl_entry);

				$index = firstidx{$_ eq $ahl_entry_field[0]} @sequenom_header;
			}


		# print "$sequenom_genotypes[0]\tgroup$ahl_entry_field[1]\t$sequenom_genotypes[$index]\n";
		print "$sequenom_genotypes[0]\tgroup$ahl_entry_field[1]\t$sequenom_genotypes[$index]\n";



		}


	}


}
