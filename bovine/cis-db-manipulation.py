#!/usr/bin/python

from __future__ import print_function
import argparse
import csv
import vcf
from Bio.Seq import Seq

## Use the csv module, as csvs are more complicated than just having a comma to delineate fields (sometimes there is a comma within a quoted field, which really borks things up)

parser = argparse.ArgumentParser(description = "Work through the CSV output of CIS-BP to adjust the coordinates to genomic coordinates", epilog = "RS Fraser, 2017-12-05")

parser.add_argument("csv", help = "The output of CIS-BP")
parser.add_argument("vcf_file", help = "A tabix indexed and bgzip compressed VCF file is required. The script will output a new field that contains the modified transcription factor binding site sequence, along with the rsid of the variant.")

args = parser.parse_args()

variant_file = vcf.Reader(open(args.vcf_file, 'r'))

with open(args.csv,'r') as file:
    fileread = csv.reader(file, delimiter = ',', quotechar = '\"')
    for entry in fileread:

        #Obtain the genomic coords from the fasta entry (1st line)
        if len(entry) == 1:
            chrom = entry[0].split(":")[2]
            chromo = "chr" + chrom
            gen_start = int(entry[0].split(":")[3])
            gen_stop = int(entry[0].split(":")[4])
            strand = int(entry[0].split(":")[5])
            gene_name = entry[0].split(":")[6]

        elif entry.count("TF ID") == 1:
            print("#chrom", "genomic_start", "genomic_end", end = "\t", sep = "\t")
            print(*entry, end = "\t", sep = "\t")
            print("var_pos", "rsID", "ref", "alt", "alt_seq", sep = "\t")
        else:
            seq = entry[5]
            tf_start = int(entry[6])
            tf_end = int(entry[7])
            tf_strand = entry[8]
            #Only consider cis acting TFs
            #Account for strand in the orig Fa file (which reports the sequence 5' - 3'; so the genome_start is the END of the sequence, and the genome_end is the START of the sequence)
            if strand == 1:
                if "F" in tf_strand:
                    #Subtract one from the genome start position because that position is actually position 0
                    tf_start_gen = (gen_start - 1) + tf_start
                    tf_end_gen = (gen_start - 1) + tf_end

                    #the use of fetch requires the vcf file to be bgziped and tabix indexed.
                    for record in variant_file.fetch(chromo, tf_start_gen,tf_end_gen):
                        newSeq = str()
                        alt_pos = tf_end_gen - record.POS + 1
                        
                        #slice and dice: in order to replace the ref base in the immutable seq variable, we need to slice up the seq var and reconstruct it.

                        #deal with the special case where alt_pos = 0 ('cause then the whole seq string gets repeated')
                        if alt_pos == 1:
                            newSeq = seq[:-alt_pos] + str(record.ALT[0]) 
                        else:
                            newSeq = seq[:-alt_pos] + str(record.ALT[0]) + seq[-alt_pos + 1:]

                        # print(alt_pos, seq, newSeq)
                        print(chromo, tf_start_gen, tf_end_gen, sep = "\t", end = "\t")
                        print(*entry, sep = "\t", end = "\t")
                        print(record.POS, record.ID, record.REF, record.ALT[0], newSeq, sep = "\t")

            elif strand == -1:
                if "F" in tf_strand:
                    #Add one from the genome start position because that position is actually position 0
                    tf_start_gen = (gen_stop + 1) - tf_end
                    tf_end_gen = (gen_stop + 1) - tf_start

                    #the use of fetch requires the vcf file to be bgziped and tabix indexed.
                    for record in variant_file.fetch(chromo, tf_start_gen,tf_end_gen):

                        newSeq = str()

                        alt_pos = record.POS - tf_start_gen + 1                   

                        #slice and dice: in order to replace the ref base in the immutable seq string variable, we need to slice up the seq var and reconstruct it.

                        #alt allele is reported from the forward strand - need to complement it. The below is a kind of ugly one liner that should probably be split up. Seq package from biopython converts the record.ALT into a seq class, which can then be complemented, and then the index 0 pulls the complemented base back out as a string.
                        compAlt = Seq(str(record.ALT[0])).complement()[0]
                        #deal with the special case where alt_pos = 0 ('cause then the whole seq string gets repeated')
                        if alt_pos == 1:
                            newSeq = seq[:-alt_pos] + compAlt 
                        else:
                            newSeq = seq[:-alt_pos] + compAlt + seq[-alt_pos + 1:]

                        # print(alt_pos, seq, newSeq)
                        print(chromo, tf_start_gen, tf_end_gen, sep = "\t", end = "\t")
                        print(*entry, sep = "\t", end = "\t")
                        print(record.POS, record.ID, record.REF, record.ALT[0], newSeq, sep = "\t")