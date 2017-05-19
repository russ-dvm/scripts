#!/usr/bin/python

from __future__ import print_function
import sys

file = sys.argv[1]
file = open(file, 'r')

for line in file:
	line = line.rstrip()
	if ">" in line:
		line_entries = line.split("|")
		for i in range(0,len(line_entries)):
			if i == 0:
				print(	line_entries[i], end="|")
			else:
				print(line_entries[i], end="\t")
		print()
