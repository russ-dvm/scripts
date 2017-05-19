#!/usr/bin/python

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(description = "This will assign the gene name a VCF file containing SNPs. Very similar to the other assign_gene name script...", epilog = "RS Fraser, 2017/04/06")

parser.add_argument("targets", help = "VCF file with the SNPs of interest")
parser.add_argument("coord", help = "File containing reference coordinates for the genes of interest. Format: gene name, genomic coordinates of the 3kb downstream that was queried.")

args = parser.parse_args()

target = open(args.targets, 'r')
coords = open(args.coord, 'r')

#For reasons unknown to me (but probably obvious) looping over the coordinate file in the "target" loop didn't really work... So instead assigned the values to a dictionary and looped over that. Whatever works...? Discovered later that probably just needed to reset the internal loop. Not worth fixing since this works.

#Create the dictionary
coord_table = {}

for gene in coords:
	if "Transcript_ID" in gene:
		continue
	else:
		gene = gene.rstrip()
		gene_fields = gene.split("\t")

		if gene_fields[5] is "1":
			coord_table[gene_fields[3]] = (gene_fields[1], gene_fields[2], gene_fields[5])
		else: 
			coord_table[gene_fields[3]] = (gene_fields[2], gene_fields[1], gene_fields[5])


#Loop over the miRNA targets. Isolate the start position of the miRNA and then cross-ref to the dictionary to identify gene name.
#Worth knowing the distance of the miRNA start position from the stop codon. Closer to stop codon may indicate that it's more likely to be actually in the 3' UTR? Not 100% sure.

for line in target:
	line = line.rstrip()
	if "##" in line:
		continue
	elif "#CHROM" in line:
		print(line, "GENE", "DISTANCE", sep="\t")
	else:
		line_fields = line.split("\t")
		target_coord = line_fields[1]
		print(line, end="\t")

		for k,v in coord_table.iteritems():
			if v[0] <= target_coord <= v[1]:
				if v[2] is "1":
					distance = int(target_coord) - int(v[0])
				else:
					distance = int(v[1]) - int(target_coord)
				if distance <= 500:
					distance_string = "<= 500"
				elif 500 < distance <= 1000:
					distance_string = "501 - 1000"
				elif 1000 < distance <= 1500:
					distance_string = "1001 - 1500"
				elif 1500 < distance <= 2000:
					distance_string = "1501 - 2000"
				elif 2000 < distance <= 2500:
					distance_string = "2001 - 2500"
				elif 2500 < distance:
					distance_string = "2501 - 3000"
				print(k, distance_string, sep="\t")
