#!/usr/bin/perl
use strict;
use warnings;
use bignum;
use diagnostics;
no warnings 'recursion';

#RSF, rfrase03@uoguelph.ca
#2016/01/06
#This script will go through the output of fisher-test-continued.pl and calculate fisher's exact test.
#THIS SCRIPT IS SUPER SLOW. LIKE REALLY SLOW. IT ALSO OUTPUTS A ONE-TAILED FISHER EXACT TEST WHICH IS NOT REALLY WHAT WE WANT... SO, I DO NOT RECOMMEND USING THIS.
#AN R-SCRIPT THAT IS MUCH FASTER AND WHICH USES A TWO-TAILED FET HAS BEEN WRITTEN AND SHOULD BE USED.
#SEE fet.R

my $usage = "perl script.pl vcf_file out_file\n";

my $vcf_file = shift || die $usage;
my $out_file = shift || die $usage;

print "YOU SHOULD REALLY USE THE R-SCRIPT INSTEAD!!!\n";
print "This might take awhile...\n";


open IN, "$vcf_file" || die "could not open $vcf_file for reading\n";
open OUT, ">$out_file" || die $usage;

print OUT "chrom\tposition\trsid\tref\talt\tp\t-logp\n";


while (<IN>) {
	s/\cM//;
	chomp;

	my $line = $_;

	if ($line =~ /chrom/) {
						#print intro lines
	}
	
	else {
		my @fields = split("\t", $line);   	#split each entry into an array			
		for (my $i = 0; $i < 5; $i++) {						
			print OUT "$fields[$i]\t";		#print first 5 fields, i.e. standard VCF fields
		}
		my $a = $fields[9];
		my $b = $fields[10];
		my $c = $fields[11];
		my $d = $fields[12];
		my $a_fac = fac($a);
		my $b_fac = fac($b);
		my $c_fac = fac($c);
		my $d_fac = fac($d);
	
		my $n = $a+$b+$c+$d;	#total count in the pop
		my $n_fac = fac($n);
	
		my $ab = $a+$b; 	#row1, i.e. ref
		my $cd = $c+$d;		#row2, i.e. alt
		my $ac = $a+$c;		#col1, i.e. normal
		my $bd = $b+$d;		#col2, i.e. diseased
	
		my $ab_fac = fac($ab);
		my $cd_fac = fac($cd);
		my $ac_fac = fac($ac);
		my $bd_fac = fac($bd);
	
		my $fet_num = $ab_fac*$cd_fac*$ac_fac*$bd_fac;
		my $fet_denom = $a_fac*$b_fac*$c_fac*$d_fac*$n_fac;
		my $fet = $fet_num/$fet_denom;
	
		my $logp = -log($fet);
	
		print OUT "$fet\t$logp\n";
# 		print "FET is $fet\n";
# 		print "log is $logp\n";
	}

	


	sub fac			#recursive factorial subroutine stolen from stack or somewhere...
	{
    my ($m) = @_;
    return 1 if($m <=1 );
    return Math::BigInt->new($m*fac($m-1));
	}
	
# 	my $n_ref_bin = 0;
# 	my $n_alt_bin = 0;
# 	my $d_ref_bin = 0;
# 	my $d_alt_bin = 0;
# 
# 		
# 	if ($line =~ /\#/) {
# 		print OUT "$line\n";					#print intro lines
# 		}
# 		
# 	else {
# 		my @fields = split("\t", $line);   	#split each entry into an array			
# 		for (my $i = 0; $i < 9; $i++) {						
# 			print OUT "$fields[$i]\t";		#print first 9 fields, i.e. standard VCF fields
# 		}
# 		my $count = $fields[9];		#variable count is assigned the value of the string of read counts
# 		my @depth = split('\"', $count);				
# 		for (my $i=0; $i<21; $i++) {
# 			if ($i==0){		#NORMAL, group1
# 				my @ref_count = split('\,', $depth[$i]);
# 				$n_ref_bin = $n_ref_bin + $ref_count[0];
# 				$n_alt_bin = $n_alt_bin + $ref_count[1];
# 				}
# 			elsif ($i>0 && $i<10) {		#DISEASED, groups11-18
# 				my @ref_count = split('\,', $depth[$i]);
# # 				print "$ref_count[0]\t";
# # 				print "$ref_count[1]\n";
# 				$d_ref_bin = $d_ref_bin + $ref_count[0];
# 				$d_alt_bin = $d_alt_bin + $ref_count[1];
# 				}
# 			elsif ($i>9 && $i<12) {			#NORMAL, group19, 2
# 				my @ref_count = split('\,', $depth[$i]);
# # 				print "$ref_count[0]\t";
# # 				print "$ref_count[1]\n";
# 				$n_ref_bin = $n_ref_bin + $ref_count[0];
# 				$n_alt_bin = $n_alt_bin + $ref_count[1];
# 				}
# 			elsif ($i>12 && $i <18) {			#NORMAL, groups3-6, 24
# 				my @ref_count = split('\,', $depth[$i]);
# # 				print "$ref_count[0]\t";
# # 				print "$ref_count[1]\n";
# 				$n_ref_bin = $n_ref_bin + $ref_count[0];
# 				$n_alt_bin = $n_alt_bin + $ref_count[1];
# 				}
# 			else {					#DISEASED, remainder
# 				my @ref_count = split('\,', $depth[$i]);
# # 				print "$ref_count[0]\t";
# # 				print "$ref_count[1]\n";
# 				$d_ref_bin = $d_ref_bin + $ref_count[0];
# 				$d_alt_bin = $d_alt_bin + $ref_count[1];
# 				}
# 			}
# 		
# 		print "$count\n";
# 		print "$n_ref_bin\t$n_alt_bin\n";
# 		print "$d_ref_bin\t$d_alt_bin\n";
# 		print OUT "$n_ref_bin\t$n_alt_bin\t$d_ref_bin\t$d_alt_bin\n";
# # 		print OUT "\n";
# 
# 	}
}		

close IN;
close OUT; 