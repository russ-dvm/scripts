#!/usr/bin/perl

use strict; use warnings;

#R. Fraser, July 14/2015
#Script to ensure that the same number of genotype entries are present on every line.
#Usage: perl <script> <input> <output>

open (my $vcf, "<", $ARGV[0]);
open (my $output, ">", $ARGV[1]);
my $num_animals = $ARGV[2];

die "Either input, output, or # of animals is missing.\nUsage: perl <correct_ploidy_var.pl> <input> <output> <num_animals>\n" if @ARGV < 3; 
die "Too many arguments.\nUsage: perl <correct_ploidy_var.pl> <input> <output> <num_animals>\n" if @ARGV > 3; 

while (my $line = <$vcf>) {
	chomp $line;
	
		
	if ($line =~ /\#/) {
		print $output "$line\n";					#print intro lines
		}
		
	else {
		
		my @fields = split("\t", $line);   				#split each entry into an array
		my @header = @fields[0 .. 8];					#create an array of the first nine columns in each line.
		shift @fields for (0 .. 8);						#alter the original array so that it only contains genotype information.
		my $genotype_num = @fields;						#total number of fields in the array now = numer of genotypes (i.e. animals).
		my @sorted = sort @fields;						#required later to handle entries in which there are more genotypes than animals. All ./. entries end up first.
		@sorted = reverse @sorted;
				
		for (my $i = 0; $i < 9; $i++) {						
			print $output "$header[$i]\t";						#print first 9 fields, i.e. standard VCF fields
		}
		
		for (my $a = 0; $a < $num_animals; $a++) {
			if ($genotype_num <= $num_animals) {				# are there less entries in the VCF files than there are animals
					if ($a < $genotype_num) {						#first, print out any existing information (this will take care of all entries that have the appropriate number)
						print $output "$fields[$a]\t";
					}
					if ($a >= $genotype_num) {						#then, enter in ./. for every missing entry until there are (9 + num of animals) total
						print $output "./.\t";
					}
				}

				elsif ($genotype_num > $num_animals) {				
						print $output "$sorted[$a]\t";				#takes the sorted array (with all ./. occurring first) and prints fields starting from the end, thus omitting extraneous ./.			
				}	
			
		}
	print $output "\n";										#add a newline to the end of each entry
	}
}	

close $vcf;
close $output;