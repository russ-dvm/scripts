#!/usr/bin/python

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(description = "Take the output of <depth_of_coverage_by_region.py> and incporate information on variation from a VCF file. Will simply add a boolean to a new column 'is_variant'.", epilog = "RS Fraser, 2017/06/14")

parser.add_argument("depth_file", help = "The output from <depth_of_coverage_by_region.py>")
parser.add_argument("vcf_file", help = "A VCF file containing variants to be annotated.")

args = parser.parse_args()


depth = open(args.depth_file, 'r')
vcf = open(args.vcf_file, 'r')

#Create list of VCF positions
vcf_list = []

for entry in vcf:
	entry = entry.rstrip()
	if "##" in entry:
		continue
	else:
		entry_fields = entry.split("\t")

		vcf_list.append(entry_fields[1])

#Parse the annotated file and see if positions match

for line in depth:
	line = line.rstrip()
	line_fields = line.split("\t")

	position = line_fields[1]

	if "chrom" in line:
		print(line, "is_variant", sep="\t")
	else:
		if position in vcf_list:
			print(line, "TRUE", sep="\t")
		else:
			print(line, "FALSE", sep="\t")