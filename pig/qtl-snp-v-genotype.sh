#!/bin/bash


#input is SNP-id


touch "$1".txt

head -1 SNP-march28.txt >> "$1".txt

grep "$1" SNP-march28.txt >> "$1".txt
grep "$1" working_gene_expression.txt >> "$1".txt