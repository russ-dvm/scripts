#!/usr/bin/perl

#Go through a sync file and remove zeros.
#Totally imperfect way of doign this, still needs a bunch of concious manip to make sure things are good.

use strict;
use warnings;


open (my $input, "<", $ARGV[0]);

while (my $line = <$input>) {
  chomp $line;

  my @full_line = split("\t", $line);

  print "$full_line[0]\t$full_line[1]\t$full_line[2]\t";

  my @pop1 = split(":", $full_line[3]);
  my @pop2 = split(":", $full_line[4]);
  my $counter_n=1;
  foreach my $value (@pop1) {
    if ($value != 0) {
      if ($counter_n == 1) {print "A:"}
      elsif ($counter_n == 2) {print "T:"}
      elsif ($counter_n == 3) {print "C:"}
      elsif ($counter_n == 4) {print "G:"}
      print "$value;";
      $counter_n = $counter_n + 1
    }
    else {$counter_n=$counter_n+1}
  }

  print "\t";

  my $counter_d = 1;

  foreach my $pop2value (@pop2) {
    if ($pop2value != 0) {
      if ($counter_d == 1) {print "A:"}
      elsif ($counter_d == 2) {print "T:"}
      elsif ($counter_d == 3) {print "C:"}
      elsif ($counter_d == 4) {print "G:"}
      print "$pop2value;";
      $counter_d = $counter_d + 1
    }
    else {$counter_d = $counter_d + 1}
  }
  print "\n";
}
