#!/usr/bin/python

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser(description = "Condenses the output of the intersection between MEME and CIS-BP by collpasing positive hits by TF Family (ie if multiple TFs of a single family bind the identical sequence they will be comma-separated into a single entry). Individual TFs will now appear as a list.", epilog = "RS Fraser, 2017-12-15")

parser.add_argument("input", help = "The intersection file.")


args = parser.parse_args()

uniqIdentifierPrev = "fake"
nameList = []
cisIDList = []
tfEnsList = []


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

        #Deal with the first line
        if "fake" in uniqIdentifierPrev: 
            nameList.append(tfName)
            cisIDList.append(cisID)
            tfEnsList.append(tfEns)


        elif uniqIdentifier in uniqIdentifierPrev:
            nameList.append(tfName)
            cisIDList.append(cisID)
            tfEnsList.append(tfEns)

        #unique identifiers don't match, meaning a new tf family is in place, so print out the old stuff
        else:
            print(*prevLine[0:7], sep = "\t", end = "\t")
            print(','.join(cisIDList), end = "\t")            
            print(','.join(nameList), end = "\t")
            print(prevLine[9], end = "\t")
            print(','.join(tfEnsList), end = "\t")
            print(*prevLine[11:], sep = "\t")
            # print(cisIDList, nameList, tfEnsList)

            #reset the lists
            nameList = []
            cisIDList = []
            tfEnsList = []

            #append current line info
            nameList.append(tfName)
            cisIDList.append(cisID)
            tfEnsList.append(tfEns)

        prevLine = lineFields
        uniqIdentifierPrev = uniqIdentifier

#Print the last line
print(*prevLine[0:7], sep = "\t", end = "\t")
print(','.join(cisIDList), end = "\t")
print(','.join(nameList), end = "\t")
print(prevLine[9], end = "\t")
print(','.join(tfEnsList), end = "\t")
print(*prevLine[11:], sep = "\t")
