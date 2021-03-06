#!/usr/bin/python

from __future__ import print_function
from collections import defaultdict
import argparse


parser = argparse.ArgumentParser(description = "Take the output of <depth_of_coverage_by_region.py> and incporate information on variation from a VCF file. Will simply add a boolean to a new column 'is_variant'.", epilog = "RS Fraser, 2017/06/14")

parser.add_argument("depth_file", help = "The output from <depth_of_coverage_by_region.py>")
parser.add_argument("vcf_file", help = "A VCF file containing variants to be annotated.")
parser.add_argument("--snpeff", help = "An optional to denote whether the vcf file has been run through SnpEff, and whether synnymous/non-synonymous info should be output in the annotation", action="store_true")
parser.add_argument("--all", help = "When set, this will output all positions in the file, and not just the variants", action = "store_true")


args = parser.parse_args()


depth = open(args.depth_file, 'r')
vcf = open(args.vcf_file, 'r')

#Create list of VCF positions
vcf_list = []

#Empty dictionary to associate a snpeff annotation to the variant position
snpeff_dict = defaultdict(list)

for entry in vcf:
	entry = entry.rstrip()
	if "#" in entry:
		continue
	else:
		entry_fields = entry.split("\t")

		# In bovine expt, there are overlapping coords on diff chroms (this wasn't a problem with equine). Needs a mech to ensure that position 15 on chr10 doesn't match position 15 on chr20. Easy solution - concatenate the chrom and position strings to make a unique identifier - eg chr10:15. This style also matches Ensembl/UCSC so I think we'll use it.

		chrom = entry_fields[0]
		var_pos = entry_fields[1]
		var_entry = chrom + ":" + var_pos
		vcf_list.append(var_entry)
		annotation_field = entry_fields[7]
		annotation_sub_field = annotation_field.split(";")
		annotation_sub_field_1 = annotation_sub_field[3]
		annotation_sub_field_2 = annotation_sub_field_1.split("|")
		annotation_sub_field_3 = annotation_sub_field_2[1]
		snpeff_dict[var_entry].append(annotation_sub_field_3)

# Parse the annotated file and see if positions match

for line in depth:
	line = line.rstrip()
	line_fields = line.split("\t")

	chrom = line_fields[0]
	position = line_fields[1]

	depth_entry = chrom + ":" + position

	if args.snpeff:
		if "chrom" in line:
			print(line, "is_variant", "variant_type", sep="\t")
		else:
			if depth_entry in snpeff_dict:
				print(line, "TRUE", snpeff_dict[depth_entry], sep="\t")
			else:
				if args.all:
					print(line, "FALSE", "FALSE", sep = "\t")
				
	else:
		if "chrom" in line:
			print(line, "is_variant", sep="\t")
		else:
			if depth_entry in snpeff_dict:
				print(line, "TRUE", sep="\t")
			else:
				if args.all:
					print(line, "FALSE", sep="\t")

