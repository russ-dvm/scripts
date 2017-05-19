#!/usr/bin/python

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(description = "Assign the name of the gene that the variant is involved in, as well as place a description group/condition in the file (optional).", epilog = "RS Fraser, 2017/05/03")

parser.add_argument("input", help = "The output from PLINK2 association test.")
parser.add_argument("-ref", nargs = "?", const = "empty", dest = "ref", help = "File containing the rsID and related gene")
parser.add_argument("-desc", nargs = '?', const = "empty", dest = 'desc_term', type=str, help= "Add a description of the subgroup being tested (e.g. enteritis, Strep zoo). Optional.")

args = parser.parse_args()

assoc = open(args.input, 'r')
ref = open(args.ref, 'r')
descriptor = args.desc_term

ref_dict = {}

for line in ref:
	line = line.rstrip()
	line_fields = line.split("\t")

	ref_dict[line_fields[2]] = line_fields[0]

for line in assoc:
	assoc_line = line.rstrip()
	if "SNP" in line and descriptor is None:
		print(assoc_line, "GENE", sep="   ")
		continue
	elif "SNP" in line and descriptor is not None:
		print(assoc_line, "GENE", "condition", sep = "    ")
		continue
	assoc_line_fields = assoc_line.split()
	gene = ref_dict[assoc_line_fields[1]]

	if descriptor is None:
		print(assoc_line, gene)
	
	else:
		print(assoc_line, gene, descriptor)