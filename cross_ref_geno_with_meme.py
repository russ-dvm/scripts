#!/usr/bin/python

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(description = "Takes a BED file with MEME coordinates of conserved upstream regions and cross-references it with an annotated version of Genomatix's SNPINSPECTOR output (tsv).", epilog = "RS Fraser, May 26, 2017")

parser.add_argument("bed", help = "BED file containing output from MEME")
parser.add_argument("tsv", help = "File from SnpInspector, annotated with genomic coordinates of the SNP.")

args = parser.parse_args()


bed_file = open(args.bed, 'r')
tsv_file = open(args.tsv, 'r')


#Read in BED as a list
bed_list = []

for line in bed_file:
	line = line.rstrip()
	line_fields = line.split("\t")
	entry_tuple = (line_fields[1], line_fields[2])
	bed_list.append(entry_tuple)


for entry in tsv_file:
	entry = entry.rstrip()
	entry_fields = entry.split("\t")
	position = entry_fields[13]
	for item in bed_list:
		if item[0] <= position <= item[1]:
			print(entry)
