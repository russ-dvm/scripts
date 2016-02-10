#!/usr/bin/perl
use strict;
use warnings;

#Ann Meyer, June 30, 2015
#This script will go through a vcf file
#and split the genotype field into a diploid field
#Modified by RSF on Sept 24 - commented out portion of the script that retained the FORMAT info.
#note that the commented out portion delivers improper formating is a genotype of 
#././././././.:0,0 is present.

my $usage = "perl script.pl vcf_file out_file\n";

my $vcf_file = shift || die $usage;
my $out_file = shift || die $usage;

open IN, "$vcf_file" || die "could not open $vcf_file for reading\n";
open OUT, ">$out_file" || die $usage;
while (<IN>) {
	s/\cM//;
	chomp;

	my $line = $_;

	if ($line =~ m/^#/) {
		print OUT "$line";
	} else {
		my @vcf_line = split('\t', $line);
		foreach my $vcf_field (@vcf_line) {
			if ($vcf_field =~ /\//) {
				my @vcf_sub = split(':', $vcf_field);
				my $field_index = 1;	
				foreach my $sub_field (@vcf_sub) {
					if ($sub_field =~ /\// && $sub_field !~ /\./) {
						my @genos = split('/', $sub_field);
						for (my $i=0; $i < scalar(@genos); $i+=2) {
							my $second = $i+1;
							if ($second == (scalar(@genos)-1) && $genos[$second] !~ /\./) {
								print OUT "$genos[$i]"."/"."$genos[$second]"."\t";
							} elsif ($genos[$second] !~ /\./) {
								print OUT "$genos[$i]"."/"."$genos[$second]\t";
							} else {
								print OUT "$genos[$i]"."/"."$genos[$second]\t";
							}
						}
						$field_index++;
					} elsif ($sub_field =~ /\// && $sub_field =~ /\./) {
						my $dots = $sub_field;
						$dots =~ s/\///g;
						for (my $i=1; $i < length($dots); $i+=2) {
							my $second = $i+1;
							unless ($second == length($dots)) {
								print OUT "./.\t";
							} else {
								print OUT "./.\t";
							}
						}						
					} #else {
						#if ($field_index == scalar(@vcf_sub)) {
						#	print OUT "$sub_field\t";
						#} else {
					#		print OUT "$sub_field".":";
					#	}
						$field_index++;
					#}
				}
			} else {
				print OUT "$vcf_field\t";
			}
		}
	}
	print OUT "\n";
}
close IN;
close OUT; 

