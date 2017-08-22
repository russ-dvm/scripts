#!/usr/bin/python

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(description = "Assign 'full' or 'seed' status to SNPs depending on where they target MREs", epilog = "RS Fraser, 2017-08-22")

parser.add_argument("seed", help = "VCF file containing seed SNPs")
parser.add_argument("full", help = "VCF file containing full-region SNPs")

args = parser.parse_args()

seed = open(args.seed, "r")
full = open(args.full, "r")

##Need both a mir list and a position list, as sometimes the same SNP affects a seed and a full
seed_list = []
mir_list = []

for line in seed:
    line = line.rstrip()
    line_fields = line.split("\t")

    seed_loc = line_fields[1]
    mir_name = line_fields[30]
    seed_list.append(seed_loc)
    mir_list.append(mir_name)


for line in full:
    line = line.rstrip()
    line_fields = line.split("\t")

    full_loc = line_fields[1]
    mir_name = line_fields[30]

    if full_loc in seed_list and mir_name in mir_list:
        print(line, "seed", sep = "\t")
    else:
        print(line, "full", sep = "\t")
