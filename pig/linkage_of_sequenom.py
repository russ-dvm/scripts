#!/bin/python

from __future__ import print_function
import argparse, sys

parser = argparse.ArgumentParser(description = "Description of the program to be written", epilog = "Name, date")

parser.add_argument("ped", help = "File containing the genotype (AA format) of the animals along with the rsIDs in a header line")
parser.add_argument("ids", help = "List of rsIDs to keep.")

args = parser.parse_args()


ped = open(args.ped, 'r')
ped_lines = ped.readlines()
sequenomedIds = ped_lines[0].rstrip().split("\t")


queryIds = open(args.ids, 'r')
queryIdList = []
outputIdList = []

for entry in queryIds:
    entry_id = entry.rstrip().split("\t")[0]
    queryIdList.append(entry_id)
    outputIdList.append(entry)

index_list = []

for query_entry in queryIdList:
    try:
        index_list.append(sequenomedIds.index(query_entry))
    except:
        pass

index_list_dedup = []
for entry in index_list:
    if entry not in index_list_dedup:
        index_list_dedup.append(entry)

outName = sys.argv[2]
outName = outName.split(".")[0]
outNameInfo = outName + ".seq.info"
outNameInfo = open(outNameInfo, 'wb')


for line in ped_lines:
    if "id" in line:
        continue
    line = line.rstrip()
    line_fields = line.split("\t")
    if line_fields[0] is "0":
        print("pig0", "pig0", 0, 0, 0, 0, sep = "\t", end = "\t")
    else:
        print(line_fields[0].replace(" ", "-"), line_fields[0].replace(" ", "-"), 0, 0, 0, 0, sep = "\t", end = "\t")
        
    for id_index in index_list_dedup:
        if "0" in line_fields[id_index]:
            print("0 0", end = "\t")
        else:
            print(line_fields[id_index][0], line_fields[id_index][1], end = "\t")
    print("")


for id_index in index_list_dedup:
    out_id = sequenomedIds[id_index]
    for full in outputIdList:
        if full.startswith(out_id):
            out = full
    outNameInfo.write(out)
