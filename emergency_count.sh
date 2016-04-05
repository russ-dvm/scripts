#!/bin/bash


for x in colec10 colec11 colec12 fcn1 fcn1-like fcn3 masp1 masp2 mbl1 mbl2 sftpa sftpd
do
	for y in 50kb_up upstream coding intron downstream
	do
	
		count=$(grep -vc \# "$y"/"$x"/*all*.vcf)
		echo -ne "$x""\t""$y""\t""$count""\n" >> emerg.txt
	done
	
done