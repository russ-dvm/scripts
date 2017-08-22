#!/usr/bin/python

from __future__ import print_function
import argparse

!parser = argparse.ArgumentParser(description = "Assign "full" or "seed" status to SNPs depending on where they target MREs", epilog = "RS Fraser, 2017-08-22")

parser.add_argument("seed", help = "VCF file containing seed SNPs")
parser.add_argument("full", help = "VCF file containing full-region SNPs")

args = parser.parse_args()

seed = open(args.seed, "r")
full = open(args.full, "r")

seed_list = []

for line in seed:
	line = line.rstrip()
	line_fields = line.split("\t")

	seed_loc = line_fields[1]
	seed_list.append(seed_loc)

print(seed_list)