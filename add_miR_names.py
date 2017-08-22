#!/usr/bin/python

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(description="Takes a VCF file and adds the name of a miR whose seed is affected by the SNP", epilog = "RS Fraser, May 2017")

parser.add_argument("vcf_file", help="A vcf input file, or at least a file with the SNP position in the 2nd column of a tab separated file")

parser.add_argument("mir_file", help="File containing miR name and target location (output from TargetScan or miRanda with genomic coords)")
args = parser.parse_args()

vcf = open(args.vcf_file, 'r')
mir = open(args.mir_file, 'r')

mir_table = {}
i = 0

for mir_line in mir:
	mir_line = mir_line.rstrip()
	mir_line_fields = mir_line.split("\t")
	mir_name_prelim = mir_line_fields[0]
	mir_name = mir_name_prelim.replace(">","")
	mir_name = mir_name + "_" + str(i)
	loc_info = mir_line_fields[1].split(":")
	strand = loc_info[5]

	#create dictionary of miRname(key) with genomic coordinates as value - note that there can be multiple hits on the same miR, so dictionary needs to be able to accomodate that. Sloppy method: add the line number to the end of each miR, thus uniquifying it. Remove the line number before printing later on (line 48).
	mir_table[mir_name] = (mir_line_fields[10], mir_line_fields[11])

	i = i + 1

for vcf_line in vcf:
	vcf_line = vcf_line.rstrip()
	if "##" in vcf_line:
		print(vcf_line)
	elif "#CHROM" in vcf_line:
		print(vcf_line, "miR", "miR_count", sep="\t")

	else:
		vcf_line_fields = vcf_line.split("\t")
		snp_pos = vcf_line_fields[1]
		mir_list = []
		for k,v in mir_table.iteritems():
			if v[0] <= snp_pos <= v[1]:
				#need to develop a list here, since a single SNP can hit more than one miR.
				fix_mir = k.split("_")
				mir = fix_mir[0]
				mir_list.append(mir)
		print(vcf_line, end = "\t")
		print(*mir_list, sep=",", end="\t")
		print(len(mir_list))
