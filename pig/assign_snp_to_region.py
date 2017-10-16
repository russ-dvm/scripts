#!/usr/bin/python

from __future__ import print_function
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(
    description="Add intervals to a VCF file",
    epilog="RS Fraser, 2017-09-28")

parser.add_argument("vcf",
    help="VCF file to annotate")

parser.add_argument("intervals",
    help="Interval file (chr:pos1-pos2)")


args = parser.parse_args()

vcf = open(args.vcf, "r")
intervals = open(args.intervals, "r")

interval_dict = defaultdict(list)

for line in intervals:
    line = line.rstrip()
    chrom, pos = line.split(":")
    pos1, pos2 = pos.split("-")
    positions = [pos1, pos2, line]
    interval_dict[chrom].append(positions)

for line in vcf:
    if "##" in line:
        continue
    elif "#" in line:
        print(line.rstrip(), "interval", sep="\t")
    else:
        line = line.rstrip()
        line_fields = line.split("\t")
        vcfChrom = line_fields[0]
        vcfPos = line_fields[1]
        possiblePos = interval_dict[vcfChrom]
        for i in possiblePos:
            if i[0] <= vcfPos <= i[1]:
                print(line, i[2], sep="\t")
                