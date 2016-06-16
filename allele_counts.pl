#!/usr/bin/perl
use strict;
use warnings;

#Ann Meyer, May 4, 2015
#This script will get the allele
#depths for the ref and alt alleles
#at each of the SNP positions
#for each of the subgroups

my $usage = "perl script.pl group_vcf out_file\n";

my $group_vcf = shift || die $usage;
my $out_file = shift || die $usage;

open IN, "$group_vcf" || die "could not open $group_vcf for reading\n";
open OUT, ">$out_file" || die $usage;

while (<IN>) {
        s/\cM//;
        chomp;

        next if (m/^#/);

        my $line = $_;

        my @snp_line = split(' ', $line);

        for (my $i=0; $i<scalar(@snp_line); $i++) {
                if ($i == 0 || $i == 1 || $i == 3 || $i == 4) {
                        print OUT "$snp_line[$i]\t";
                } elsif ($i > 8) {
                        if ($snp_line[$i] =~ /:/ && $snp_line[8] =~ /AD/) {
                                my @sub_group = split(':', $snp_line[$i]);
                                my ($ref, $alt) = split /,/, $sub_group[1];
                                print OUT "$ref\t$alt\t";
                        } elsif ($snp_line[$i] =~ /:/ && $snp_line[8] =~ /DP/) {
                                my @sub_group = split(':', $snp_line[$i]);
                                my $dp = $sub_group[1];
                                print OUT "$dp\t0\t";
                        } else {
                                print OUT "0\t0\t";
                        }
                }
        }
        print OUT "\n";
}
close IN;
close OUT;
