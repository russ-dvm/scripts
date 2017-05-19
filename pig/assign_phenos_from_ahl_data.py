#!/usr/bin/python


from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(description="Using information found in the AHL Data spreadsheet, this creates several PED files that assign phenotypes for the various pathogens or conditions that were tested for. Healthy animals are controls, animals that test positive are cases, and diseased animals that were not tested or test negative are not considered.", epilog = "RS Fraser, April 25, 2017")

parser.add_argument("ped_file", help="The original PED file that is used as a template")
parser.add_argument("ahl_file", help="The AHL file containing the pathogen/condition data. Pigs as columns, conditions as rows. Test positive = 1; anything else is set as a phenotype to ignore")

args = parser.parse_args()

ped = open(args.ped_file, 'r')
ahl = open(args.ahl_file, 'r')

x_dict = {}

for line in ahl:

	line = line.rstrip()
	if "animal_id" in line:
		pig_id_keys = line.split("\t")

	else:
		line_fields = line.split("\t")
		rsid_dict = line_fields[0] + "_dict"
		rsid_dict = {}

		for i in range(len(line_fields)):
			if i == 0:
				rsid = line_fields[i]
			else:
				rsid_dict[pig_id_keys[i]] = line_fields[i]
			
		x_dict[rsid] = rsid_dict


for k,v in x_dict.items():

	filename_str = k + ".ped"
	filename = open(filename_str, 'w')

	for ped_line in ped:

		ped_line = ped_line.rstrip()
		ped_line_fields = ped_line.split("\t")

		if "-" not in ped_line:
			print(ped_line, file = filename)
		
		else:
		
			for k2, v2 in x_dict[k].items():
		
				if ped_line_fields[0] == k2:

					for i in range(5):
						print(ped_line_fields[i], sep="\t", end="\t", file = filename)
					
					if "1" in x_dict[k][k2]:
						print("2", sep="", end="\t", file = filename)

					else:
						print("-9", end="\t", file = filename)

					for i in range(6, len(ped_line_fields)):
						print(ped_line_fields[i], sep="\t", end="\t", file = filename)

					print("", file = filename)
					
	#holy shit did it take awhile to debug this bastard. The issue with this nested for loop is that once the ped file has been iterated over, python doesn't reset it back to the start... so the first key value from teh dictionary would iterate over the entire ped file, and the ped file would be at it's "end". The seek function resets the file's position to the defined value - thus, for every key iteration, it starts the ped file back at the start. Whoosh....
	ped.seek(0)
