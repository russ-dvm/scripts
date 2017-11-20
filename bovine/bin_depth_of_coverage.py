#!/usr/bin/python

from __future__ import print_function
import argparse


parser = argparse.ArgumentParser(description = "Takes the output of Depth of Coverage (from GATK) - where the depth of each base is reported - and converts into specified bins", epilog = "Name, date")

parser.add_argument("input", help = "Output of GATK Depth of Coverage")
parser.add_argument("binSize", help = "Integer value for the desired bin size")

args = parser.parse_args()

doc = open(args.input, 'r')
bins = int(args.binSize)


sampleNum = 0
depth = 0
i = 0

for line in doc:
    line = line.rstrip()
    line_fields = line.split("\t")

    #Determine number of samples (so it's not experiment specific)
    if "Locus" in line:
        for field in line_fields:
            if "Depth_for" in field:
                sampleNum = sampleNum + 1

    else:
        chrom = line_fields[0].split(":")[0]
        pos = line_fields[0].split(":")[1]
        #Determine start of the interval
        if i == 0:
            startChrom = chrom
            startPos = pos

        #old pos won't be defined for the first line - use try to workaround
        try:
            if int(pos) == int(oldPos) + 1:

                depth = depth + float(line_fields[1])
                
                i = i + 1

                #Print the results at the end of the determined interval length
                if i == bins:
                    endPos = pos
                    avgDepth = depth/sampleNum/bins

                    print(chrom, startPos, endPos, avgDepth)

                    #Reset the variables
                    i = 0
                    depth = 0
            else:
                # depth = depth + float(line_fields[1])
                endPos = oldPos
                test = int(endPos) - int(startPos)
                avgDepth = depth/sampleNum/test
                print(oldChrom, startPos, endPos, avgDepth)
                i = 0
                depth = 0

        #oldPos should only be undefined for the first line; after that it will always be defined. so the exception is the first line.
        except:
            pass
            
        oldPos = pos
        oldChrom = chrom
