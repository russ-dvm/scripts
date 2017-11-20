#!/usr/bin/python

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(description = "Prepare a VCF file containing only SNPs of interest for analysis with Polyphen2.", epilog = "RS Fraser, 2017/04/06")

parser.add_argument("vcf", help = "VCF file with the SNPs of interest")

args = parser.parse_args()

vcf = open(args.vcf, 'r')

aa_dict = {"Arg":"R", "His":"H", "Lys":"K", "Asp":"D", "Glu":"E", "Ser":"S", "Thr":"T", "Asn":"N", "Gln":"Q", "Cys":"C", "Sec":"U", "Gly":"G", "Pro":"P", "Ala":"A", "Val":"V", "Ile":"I", "Leu":"L", "Met":"M", "Phe":"F", "Tyr":"Y", "Trp":"W"}

for line in vcf:
	if "#" in line:
		continue
	else:
		line = line.rstrip()
		if "synonymous_variant" in line:
			continue
		else:
			line_fields = line.split("\t")
			annotation_fields = line_fields[7].split("|")
			transcript = annotation_fields[6].split(".")
			transcript = transcript[0]
			protein_info = annotation_fields[10].split(".")
			orig_aa = aa_dict[str(protein_info[1][:3])]
			mutated_aa = aa_dict[str(protein_info[1][-3:])]
			pos = int(protein_info[1][3:-3])

			print(transcript, pos, orig_aa, mutated_aa, line_fields[2], line_fields[3], line_fields[4])