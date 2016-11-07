#!/usr/bin/python

from __future__ import print_function
import sys

with open(sys.argv[1]) as file:
    lines = file.readlines()

for i in range(1, len(lines)):
    line_fields = lines[i].split()

    i = i - ((int(line_fields[0])-1) * 93)

    if i % 93 == 0:
        print(line_fields[2],"double", "blank", "blank", sep = "\t", end = "\n")
        print("")
        print("")
    elif i % 8 != 0:
        print(line_fields[2], end="\t")
    elif i % 8 == 0:
        print(line_fields[2])
