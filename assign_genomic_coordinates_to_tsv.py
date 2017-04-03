#!/usr/bin/python

## RSF, 2017/03/31
## Retrieve genomic coordinates for the SNPs used in the SNP_inspector analysis of TFBS

from __future__ import print_function
import sys

#open the "*.geno.txt" file that contains the genomic coordinates of the gene
coord = sys.argv[1]
coord = open(coord, 'r')

#open the TSV file output from SNP_inspector
tsv = sys.argv[2]
tsv = open(tsv, 'r')

#Obtain the chromosome number, strand, and the genomic coordinate corresponding to the start of the region analyzed (usually about 3kb upstream of the Ensembl gene coordinates)
for coord_line in coord:
    if ">" in coord_line:
        coord_line_fields = coord_line.split("|")
        for coord_line_field in coord_line_fields:
            if "start" in coord_line_field:
                start_field = coord_line_field.split("=")
                start_pos = start_field[1]
            elif "end" in coord_line_field:
                end_field = coord_line_field.split("=")
                end_pos = end_field[1]
            elif "chr" in coord_line_field:
                chr_field = coord_line_field.split("=")
                chr_num = chr_field[1]
            elif "str" in coord_line_field:
                str_field = coord_line_field.split("=")
                if "(-)" in str_field[1]:
                    strand = "reverse"
                else:
                    strand = "forward"

#Open the TSV file and determine the "allele position" and use that to recalculate the actual genomic position of the SNP. Note that it differs if the gene is on the forward or reverse strand.
for tsv_line in tsv:
    tsv_line = tsv_line.rstrip()

    if "Seq" not in tsv_line:
        if strand is "forward":
            tsv_line_fields = tsv_line.split("\t")
            pos = int(tsv_line_fields[1]) + int(start_pos) -1
            print(tsv_line, chr_num, pos, sep="\t")
        elif strand is "reverse":
            tsv_line_fields = tsv_line.split("\t")
            pos = int(end_pos) - int(tsv_line_fields[1]) + 1
            print(tsv_line, chr_num, pos, sep="\t")
    else:
        print(tsv_line, "chrom", "position", sep="\t")
