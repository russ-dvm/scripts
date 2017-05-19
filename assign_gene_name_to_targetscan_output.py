#!/usr/bin/python

#RSF, 2017/04/06

#This will assign the gene name to the output of the TargetScan file. Requires the genomic positions of the miRNA targets to have been added with the "get_genomic_coords_from_Targetscan.py" script.

from __future__ import print_function
import sys

#Input is two files: the targetscan output file (with genomic coordinates appended) and a file with the gene name, and genomic coordinates of the 3kb downstream that was queried.

target = sys.argv[1]
target = open(target, 'r')

coords = sys.argv[2]
coords = open(coords, 'r')

#For reasons unknown to me (but probably obvious) looping over the coordinate file in the "target" loop didn't really work... So instead assigned the values to a dictionary and looped over that. Whatever works...?

#Create the dictionary
coord_table = {}

for gene in coords:
	if "Transcript_ID" in gene:
		continue
	else:
		gene = gene.rstrip()
		gene_fields = gene.split("\t")

		if gene_fields[5] is "1":
			coord_table[gene_fields[3]] = (gene_fields[1], gene_fields[2])
		else: 
			coord_table[gene_fields[3]] = (gene_fields[2], gene_fields[1])


#Loop over the miRNA targets. Isolate the start position of the miRNA and then cross-ref to the dictionary to identify gene name.
#Worth knowing the distance of the miRNA start position from the stop codon. Closer to stop codon may indicate that it's more likely to be relevant? Not 100% sure.

for line in target:
	line = line.rstrip()
	line_fields = line.split("\t")
	target_coord = line_fields[14]
	print(line, end="\t")

	for k,v in coord_table.iteritems():
		if v[0] <= target_coord <= v[1]:
			distance = int(target_coord) - int(v[0])
			if distance <= 500:
				distance_string = "<= 500"
			elif 501 < distance <= 1000:
				distance_string = "501 - 1000"
			elif 1001 < distance <= 1500:
				distance_string = "1001 - 1500"
			elif 1501 < distance <= 2000:
				distance_string = "1501 - 2000"
			elif 2001 < distance <= 2500:
				distance_string = "2001 - 2500"
			elif 2501 < distance:
				distance_string = "2501 - 3000"
			print(k, distance_string, sep="\t")
