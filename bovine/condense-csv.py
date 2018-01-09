#!/usr/bin/python

from __future__ import print_function
from collections import defaultdict
import argparse
import csv

parser = argparse.ArgumentParser(description = "Condenses the output of CIS-BP .csv file by collpasing positive hits by TF Family (ie if multiple TFs of a single family bind the identical sequence they will be comma-separated into a single entry). Individual TFs will now appear as a list.", epilog = "RS Fraser, 2018-01-05")

parser.add_argument("input", help = "The CSV file.")


args = parser.parse_args()

uniqIdentifierPrev = "fake"

preList = []
cisIDList = []
NameList = []
midList = []
EsnIDList = []
afterList = []

# infoDict = defaultdict(preList, cisIDList, NameList, midList, EsnIDList, afterList)
infoDict = defaultdict(list)


with open(args.input, 'r') as file:
    fileread = csv.reader(file, delimiter = ',', quotechar = '\"')
    for lineFields in fileread:
        if "chrom" in lineFields[0] or "TF ID" in lineFields[0]:
            continue
        else:
            family = lineFields[4]
            tfEns = lineFields[3]
            motifID = lineFields[2]
            tfName = lineFields[1]
            sequence = lineFields[5]
            start = lineFields[6]
            stop = lineFields[7]
            uniqIdentifier = family+":"+start+":"+stop+":"+sequence

            if uniqIdentifier not in infoDict:
                infoDict.setdefault(uniqIdentifier, [[],[],[],[],[],[]])
                # infoDict[uniqIdentifier][0].extend(lineFields[0:2])
                infoDict[uniqIdentifier][2].append(lineFields[2])
                infoDict[uniqIdentifier][5].extend(lineFields[4:])
            # infoDict[uniqIdentifier] = [lineFields[0:7], cisIDList,[], lineFields[9], [], lineFields[11:]]
            # infoDict[uniqIdentifier][0].append(motifID)
            infoDict[uniqIdentifier][0].append(lineFields[0])
            infoDict[uniqIdentifier][1].append(tfName)
            infoDict[uniqIdentifier][3].append(tfEns)

    for k,v in infoDict.items():
        for i in range(len(v)):
            if i == 0 or i == 1 or i == 3:
                print(','.join(v[i]), end = "\t")
            elif i == 5:
                print('\t'.join(v[i]))
            else:
                print('\t'.join(v[i]), end = "\t")