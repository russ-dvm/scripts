#!/usr/bin/python

from __future__ import print_function
import sys

coords = sys.argv[1]
coords = open(coords, 'r')

def get_position(entry):
	entry_fields = entry.split("=")
	position = entry_fields[1]
	return position

for line in coords:
	if ">" in line:
		line = line.rstrip()
		line_fields = line.split("|")
		chrom = get_position(line_fields[6])
		start_coord = get_position(line_fields[9])
		end_coord = get_position(line_fields[10])
		print("chr", end="")
		print(chrom, start_coord, end_coord, sep="\t")
