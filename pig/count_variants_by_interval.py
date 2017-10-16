#!/usr/bin/python

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(description = "Count the number of variants intersecting bins of an interval file", epilog = "RS Fraser, 2017-10-13")

parser.add_argument("interval", help = "Interval file")
parser.add_argument("vcf", help = "VCF File")

args = parser.parse_args()

interval = open(args.interval, 'r')
vcf = open(args.vcf,'r')

vcfList = []

for line in vcf:
    if "#" in line:
        continue
    chr_pos = [line.split("\t")[0], int(line.split("\t")[1])]
    vcfList.append(chr_pos)

for line in interval:
    line = line.rstrip()
    chrom = line.split(":")[0]
    start = int(line.split(":")[1].split("-")[0])
    end = int(line.split(":")[1].split("-")[1])
    counter = 0
    for entry in vcfList:
        if chrom in entry[0] and start <= entry[1] <= end:
            counter = counter + 1
    print(line, counter, sep ="\t")