#!/usr/bin/perl

use strict; use warnings;

#R. Fraser
#Used to cleanup the files generated by make_diploid into a Plink compatible input.

open (my $vcf, "<", $ARGV[0]);
open (my $output, ">", $ARGV[1]);

while (my $line = <$vcf>) {
	chomp $line;
	
		
		my @fields = split("\t", $line);   				#split each entry into an array
		my $fieldsnum = @fields;
		print "$fieldsnum\n";
}	

close $vcf;
close $output;