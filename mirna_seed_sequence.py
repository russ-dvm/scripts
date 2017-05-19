#!/usr/bin/python

from __future__ import print_function
import sys

#RSF, April 4, 2017.
#Adapt mature miRNA sequences obtained from the miRDB (on April 3, 2017) to the format required by TargetScan.

#Format of the output file should be miRNA family name, seed sequence, and taxonomic ID (in tab separated format). Equine taxonomic ID = 9796.


mirna = sys.argv[1]
mirna = open(mirna, 'r')

for mirna_entry in mirna:
    mirna_entry = mirna_entry.rstrip()
    if ">" in mirna_entry:
        mirna_entry_line = mirna_entry.split(" ")
        mirna_family_name = mirna_entry_line[0].replace(">", "")
        print(mirna_family_name, sep="\t", end="\t")
    else:
        print(mirna_entry[1:8], "9796", sep="\t")
