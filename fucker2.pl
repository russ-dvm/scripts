#!/usr/bin/perl

use strict;
use warnings;

open (my $pivot, "<", $ARGV[0]);

my @pivot_array = <$pivot>;
my (@geno_array);

print "ID	group1	group2	group3	group4	group5	group6	group7	group8	group9	group10	group11	group12	group13	group14	group15	group16	group17	group18	group19	group20	group21	group22	group23	group24	group25	group26	group27	group28	group29	group30	group31	group32	group33	group34	group35	group36	group37	group38	group39	group40	group41	group42	group43	group44	group45	group46	group47	group48	group49	group50	group51	group52	group53	group54	group55	group56	group57	group58	group59	group60	group61	group62	group63	group64	group65	group66	group67	group68	group69	group70	group71	group72	group73	group74	group75	group76	group77	group78	group79	group80	group81	group82	group83	group84	group85	group86	group87	group88	group89	group90	group91	group92	group93	group94	group95	group96	group97	group98	group99	group100	group101	group102	group103	group104\n";

for (my $i=0; $i < scalar @pivot_array; $i++) {

	chomp $pivot_array[$i];
	my @pivot_field = split(/\t/, $pivot_array[$i]);

	push @geno_array, $pivot_field[2];

	if ($i == 103 || ($i > 200 && ($i + 1) % 104 == 0)) {
		local $" = "\t";
		print "$pivot_field[0]\t@geno_array\n";
		@geno_array=();
	}
}