#!/usr/bin/perl
use strict;
use warnings;

#RSF
#Script to get chromosomal positions for the allelic association output of Plink

my $usage = "perl script.pl file_with_positions allele_file out_file\n";

open (my $vcf_with_pos, "<", $ARGV[0]);
open (my $allele_file, "<", $ARGV[1];
open (my $output, ">", $ARGV[2]);

die "Usage: perl <script> <vcf file with positions> <allele file>" unless @ARGV == 3;


my %positions;

while (my $line = <$vcf_with_pos>) {
        s/\cM//;
        chomp $line;
		my @vcf_fields = split("\t", $line);
        $positions{$vcf_fields[0]}{$vcf_fields[1]}++;

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

