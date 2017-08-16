#!/usr/bin/python

from __future__ import print_function
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description = "Annotates the output from GATKs DepthOfCoverage utility into exons, introns, upstream, and downstream. Requires a list of coordinates for the features you are interested in..", epilog = "RS Fraser, 2017/06/09")

parser.add_argument("-depth", dest="gatk_file", help = "Output file from GATK DepthOfCoverage")
parser.add_argument("-feature", dest = 'feature_file', help = "File should be organized more or less like a BED file: chr, start, stop. Subsequent columns must include a <feature> column that identifies the the region as exon, intron, upstream, etc. If fields gene_id, transcript_id, protein_id are present they will be output.")

args = parser.parse_args()

gatk = open(args.gatk_file, "r")
exons = open(args.feature_file, "r")

#read the exon info into a dictionary
##cool - defaultdict allows you to create a dictionary without default key values, very useful for downstream appending. https://stackoverflow.com/questions/3419147/appending-values-to-dictionary-in-python. Requires the import above...

exon_dict = defaultdict(list)

for entry in exons:
	entry = entry.rstrip()
	if "feature" in entry:
		header_list = entry.split("\t")
		gene_id_index = header_list.index("gene_id")
		tx_id_index = header_list.index("transcript_id")
		pr_id_index = header_list.index("pr_id")
		gene_name_index = header_list.index("gene_name")
		feature_index = header_list.index("feature")
		feature_id_index = header_list.index("feature_id")
		strand_index = header_list.index("strand")
	else:
		entry_fields = entry.split("\t")
		quick = [entry_fields[1], entry_fields[2], entry_fields[gene_id_index], entry_fields[tx_id_index], entry_fields[pr_id_index], entry_fields[strand_index], entry_fields[feature_index], entry_fields[feature_id_index], entry_fields[gene_name_index]]
		exon_dict["chr" + entry_fields[0]].append(quick)


pos_list = []

for line in gatk:
	line = line.rstrip()

 	if "Locus" in line:
 		print("chrom", "pos", "depth", "gene_name", "feature", "feature_id", "gene_id", "transcript_id", "protein_id", "strand", sep="\t")
 	else:
	 	line_fields = line.split("\t")
	  	chrom_pos = line_fields[0].split(":")
	 	chrom = chrom_pos[0]
	 	position = int(chrom_pos[1])
	 	depth = int(line_fields[1])

	 	for dict_chrom,list_of_features in exon_dict.iteritems():
	 		#does the chromosome match
	 		if chrom == dict_chrom:
	 			for feature_entry in list_of_features:
	 				if "na" in feature_entry[0] or "na" in feature_entry[1]:
	 					continue
	 				elif int(feature_entry[0]) <= position <= int(feature_entry[1]):
	 					 print(dict_chrom, position, depth, feature_entry[8], feature_entry[6], feature_entry[7], feature_entry[2], feature_entry[3], feature_entry[4], feature_entry[5], sep="\t")

	 					#pos_list.append(position)


##Want to find the values that are missing? Uncomment this along with the pos_list stuff above and have a look. This takes forever meaning there is probably a much better way to do it...

#gatk.seek(0)

#for line2 in gatk:

#	fields = line2.split("\t")
#	if "Locus" in line2:
#		continue
#	else:
#	  	chrom_pos2 = fields[0].split(":")
#	 	posn = int(chrom_pos2[1])
#
#	if posn in pos_list:
#		continue
#	else:
#		print(posn)