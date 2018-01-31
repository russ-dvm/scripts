#!/usr/bin/python

from __future__ import print_function
import sys

with open(sys.argv[1], 'r') as vcf:
    for line in vcf:
        line = line.rstrip()
        if "##fileformat" in line:
            print(line)
        if "##INFO=<ID=AF" in line:
            print(line)
        if "##contig" in line:
            print(line)
        if "##reference" in line:
            print(line)
        if "#CHROM" in line:
            head = line.split("\t")
            print(*head[0:8], sep = "\t")
        if "#" not in line:
            lineFields = line.split("\t")
            keepFields = lineFields[0:7]
            info = lineFields[7].split(";")
            print(*keepFields, sep = "\t", end = "\t")
            print(info[1])

