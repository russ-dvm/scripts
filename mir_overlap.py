#!/usr/bin/python

from __future__ import print_function
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description = "Take the output from miRanada and Targetscan and identify intersecting miRs by region.", epilog = "RS Fraser, 2017-08-21")

parser.add_argument("-mir", "--miRanada", dest = "mir_file", help = "The output from miranada", required = True)
parser.add_argument("-ts", "--targetscan", dest = 'targ_scan', help = "The output from targetscant.", required = True)

args = parser.parse_args()

mir = open(args.mir_file, 'r')
target = open(args.targ_scan, 'r')

#mir file has a leading ">", then mir name, then fasta info of the gene
#targetscan file has fasta info of the gene, then mir name,
#both are tab delimited

#quick check on the intersection of miRs between genes - tspecific coordinates irrelevant

#create a dict with key as a gene/fasta info and list of miRs as values

mir_dict = defaultdict(list)

for line in mir:
    line = line.rstrip()
    line_fields = line.split("\t")
    fasta = line_fields[1]
    mir = line_fields[0]
    mir_fields = mir.split(">")
    mir_name = mir_fields[1]
    mir_dict[fasta].append(mir_name)

for line in target:
    line_fields = line.split("\t")
    ts_fasta = line_fields[0]
    ts_mir = line_fields[1]

    for k,v in mir_dict.items:
        print(k)
