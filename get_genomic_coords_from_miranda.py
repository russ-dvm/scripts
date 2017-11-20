#!/usr/bin/python

#Translate the coordinates of the miR targets from miRanda to genomic coordinates
#The direct output of miRanda has a lot of information and is built in a complicated fashion. The relevant information is preceded by a ">", so a grep ">" directed to a new file provides a nice summary. This script is designed to deal with that summary.

#The info in the summarized files is kind of doubled (there are both > and >> lines - that contain similar but not identical info). The start/stop of the miR is contained in the ">" line. Start/stop coordinates are within field 2.

#Because this aligns the full miR to the ref, it doesn't provide coordinates for the seed sequence. This script adds those in by taking the 

from __future__ import print_function
import sys

results = sys.argv[1]
results = open(results, 'r')

for line in results:
	if ">>" in line:
		continue
	else:
		line = line.rstrip()
		line_fields = line.split("\t")
		fasta_header = line_fields[1]

		#find genomic start and end of the gene
		fasta_header_fields = fasta_header.split(":")
		strand = fasta_header_fields[5]
		chrom = fasta_header_fields[0]
		lower_genome_coord = fasta_header_fields[4]
		higher_genome_coord = fasta_header_fields[5]

		mir_coords = line_fields[5].split(" ")
		mir_start = mir_coords[0]
		mir_end = mir_coords[1]

		#return genomic coordinates, and determine 5' vs 3' (depends on strand)
		if strand is "1":
			five_coord = lower_genome_coord
			genomic_start = int(five_coord) + int(mir_start)
			genomic_end = int(five_coord) + int(mir_end)

		else:
			five_coord = higher_genome_coord
			genomic_start = int(five_coord) - int(mir_end)
			genomic_end = int(five_coord) - int(mir_start) 


		#return seed sequence coordinates
		#it seems miRanda always returns one extra base at the 3' end of the query
		#Therefore seed start should = the genomic end + 1 bp + 6 bp
		if strand is "1":
			five_coord = lower_genome_coord
			seed_start = int(five_coord) + int(mir_end) - 8
			seed_end = int(five_coord) + int(mir_end) - 2

		else:
			five_coord = higher_genome_coord
			seed_start = int(five_coord) - int(mir_end) + 2
			seed_end = int(five_coord) - int(mir_end) + 8


		print(line, chrom, genomic_start, genomic_end, seed_start, seed_end, sep="\t")
