#!/usr/bin/perl
use strict;
use warnings;

#Ann Meyer, May 12, 2015
#This script will match
#SNPs from the 50K chip
#to GATK SNPs

my $usage = "perl script.pl snpchip_file gatk_snp_file out_file\n";

my $snpchip = shift || die $usage;
my $gatk_SNPs = shift || die $usage;
my $out_file = shift || die $usage;

my %SNPs;
open IN, "$snpchip" || die "could not open $snpchip for reading\n";
while (<IN>) {
        s/\cM//;
        chomp;

        my ($chrom, $pos) = split /\t/;

        $SNPs{$chrom}{$pos}++;

}
close IN;

open OUT, ">$out_file" || die $usage;
open IN, "$gatk_SNPs" || die "could not open $gatk_SNPs for reading\n";
while (<IN>) {
        s/\cM//;
        chomp;

        my ($chrom, $pos) = split /\t/;

        if ($SNPs{$chrom}{$pos}) {
                print OUT "$chrom\t$pos\n";
        }
}
close IN;
close OUT;

