#!/usr/bin/python

from __future__ import print_function
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description = "Condenses the output of the intersection between MEME and CIS-BP by collpasing positive hits by TF Family (ie if multiple TFs of a single family bind the identical sequence they will be comma-separated into a single entry). Individual TFs will now appear as a list.", epilog = "RS Fraser, 2017-12-15")

parser.add_argument("input", help = "The intersection file.")


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
    for line in file:
        line = line.rstrip()
        lineFields = line.split("\t")

        family = lineFields[11]
        tfEns = lineFields[10]
        cisID = lineFields[7]
        tfName = lineFields[8]
        sequence = lineFields[12]
        rsID = lineFields[18]

        uniqIdentifier = family+":"+rsID+":"+sequence

        if uniqIdentifier not in infoDict:
            infoDict.setdefault(uniqIdentifier, [[],[],[],[],[],[]])
            infoDict[uniqIdentifier][0].extend(lineFields[0:7])
            infoDict[uniqIdentifier][3].append(lineFields[9])
            infoDict[uniqIdentifier][5].extend(lineFields[11:])
        # infoDict[uniqIdentifier] = [lineFields[0:7], cisIDList,[], lineFields[9], [], lineFields[11:]]
        infoDict[uniqIdentifier][1].append(cisID)
        infoDict[uniqIdentifier][2].append(tfName)
        infoDict[uniqIdentifier][4].append(tfEns)

for k,v in infoDict.items():
    for i in range(len(v)):
        if i == 1 or i == 2 or i == 4:
            print(','.join(v[i]), end = "\t")
        elif i == 5:
            print('\t'.join(v[i]))
        else:
            print('\t'.join(v[i]), end = "\t")