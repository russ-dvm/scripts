use strict;
use warnings;

#Find ungenotyped samples in the sequenom data.
#RFraser

open (my $sequenom, "<", $ARGV[0]);

my @sequenom_array = <$sequenom>;
chomp(@sequenom_array);

my @sequenom_header_fields = split(/\t/, $sequenom_array[0]);

foreach my $sequenom_entry (@sequenom_array) {

	if ($sequenom_entry =~ /family_id/){}
	else {

		my @sequenom_fields = split(/\t/, $sequenom_entry);
		print "$sequenom_fields[0]\t";

		for (my $i = 0; $i < scalar @sequenom_fields; $i++) {

			if ($i > 5 && $sequenom_fields[$i] =~ /\d/){
			print "$sequenom_header_fields[$i]\t";
			}
		}
	print "\n";

	}

}