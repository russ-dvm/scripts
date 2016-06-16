use strict;
use warnings;

#Identify which groups were not genotyped for each SNP in a vcf file.

#RSF, rfrase03@uoguelph.ca


open (my $vcf, "<", $ARGV[0]);

my @vcf_array = <$vcf>;
my (@header_array);

foreach my $entry (@vcf_array) {
	chomp $entry;
	if ($entry =~ /\##/){}
	elsif ($entry =~ /#CHROM/) {

		@header_array = split(/\t/, $entry);
	}

	else {

		my @entry_fields = split(/\t/, $entry);
		for (my $i = 0; $i < scalar @entry_fields; $i++) {

			if ($i == 2){
				print "$entry_fields[$i]\t";
			}

			if ($entry_fields[$i] =~ /\.\/\./) {
				print "$header_array[$i]\t";
			}
		}
		print "\n";
	}

}