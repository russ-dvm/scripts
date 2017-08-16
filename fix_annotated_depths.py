#!/bin/ptyhon

from __future__ import print_function
import sys



file = open(sys.argv[1], 'r')

for line in file:
	line = line.rstrip()
	line_fields = line.split("\t")
	feature = line_fields[4]
	if "downstream" in feature:
		for item in line_fields:
			print(item)
		print(line_fields[:11], "downstream", sep = "\t")
	elif "upstream" in feature:
		for item in line_fields:
			print(item, sep ="\t", end="\t")
			print("upstream", sep = "")
	elif "intron" in feature:
		print(line_fields[:11], "intron", sep = "\t")
	else:
		print(line_fields)
