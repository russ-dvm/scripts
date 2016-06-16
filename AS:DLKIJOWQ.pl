#!/usr/bin/perl

use strict;
use warnings;


open (my $temp, "<", $ARGV[0]);

my @temp_array = <$temp>;

for (my $i=0; $i < scalar @temp_array; $i++){
	chomp $temp_array[$i];

	my $temp_field = split(/\t/, $temp_array[$i]);

	if ($i % 3 = 0) {
		print "$temp_field[0]\ttemp_field[1]\t";
	}
	else {
		print "$temp_field[3]";
	}
	print "\n";

}