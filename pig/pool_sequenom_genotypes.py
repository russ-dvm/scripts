#!/usr/bin/python

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(description = "Takes the output from pig_ahl_ids_to_group_genotypes.pl and reformts it into a style acceptable for pooling_vs_sequenom. This means that each line will have an rsid, and hexaploid genotypes for all the groups.", epilog = "RS Fraser, April 27, 2017.")

parser.add_argument("input", help = "The output from pig_ahl_ids_to_group_genotypes.pl")
parser.add_argument("-ref", dest = 'ref_file', help = "A tab-separated reference file with rsIDs and the reference and alternate allele. Optional.")

args = parser.parse_args()

genos = open(args.input, 'r')
ref = open(args.ref_file, 'r')

ref_dict = {}

for ref_line in ref:
	if "rs" not in ref_line:
		continue
	else:
		ref_line = ref_line.rstrip()
		ref_line_fields = ref_line.split("\t")
		ref_dict[ref_line_fields[0]] = (ref_line_fields[1],ref_line_fields[2])

i = 1
hexaploid = ""
geno_list = []
text_list = []
group_list = []



for line in genos:
	line = line.rstrip()
	line_fields = line.split("\t")


	genotype = line_fields[2]
	hexaploid = hexaploid + genotype
	

	if i % 3 == 0:
		geno_list.append(hexaploid)
		hexaploid = ""
		group_id = line_fields[1]
		group_list.append(group_id)

	if i == 306:
		print("ID", "REF", "ALT", *group_list, sep="\t")

	if i % 306 == 0:
		print(line_fields[0], sep="\t", end="\t")
		try:
			for thing in ref_dict[line_fields[0]]:
				print(thing, end="\t")
		except KeyError:
			print("Missing", "Missing", sep="\t", end="\t")
		print(*geno_list, sep="\t")
		geno_list = []

	i = i + 1

