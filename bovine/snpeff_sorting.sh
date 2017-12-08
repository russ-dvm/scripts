#!/bin/bash

#RSF, rfrase03@uoguelph.ca
#runs snpEff and the sorts all of the SNPs into directories based on gene and region. sweeeeet.


#variables
input=$1
base=$2
interval=/media/russ/data/bovine/ref_files/intervals.txt
usage=$(echo "snpeff_sorting.sh <input> <output>")
# snpeff_genome=UMD3.1.82
snpeff_genome=UMD3.1.90.custom
gatk_genome=/media/russ/data/bovine/genome/Sequence/WholeGenomeFasta/genome.chr.fa

#check to see if the command has been issued properly - requires input and output
if [[ $input == "-usage" ]]
then
	echo
	echo $usage
	echo
elif [[ -z "$1" ]] || [[ -z "$2" ]]
then
	echo
	echo "###########################ERRROR###########################"
	echo "Incorrect syntax - one of either input or output is missing."
	echo "Correct usage:"
	echo -e "\t$usage"
	echo "############################################################"
	echo
else
	read -p "Run snpeff (y or Y)? " run_snpeff
	if [[ $run_snpeff == "y" ]] || [[ $run_snpeff == "Y" ]] || [[ $run_snpeff == "yes" ]]
	then
	##First run snpEff at defaults##
	echo
	echo "Running snpEff using $snpeff_genome, upstream/downstream 5kb."
	java -Xmx10G -jar ~/java/snpEff/snpEff.jar eff -noStats $snpeff_genome $input > $base.vcf
	##Run snpEeff and defined upstream and downstream regions as 50kb from AUG##
	echo
	echo "Running snpEff using $snpeff_genome, upstream/downstream 50kb."
	java -Xmx10G -jar ~/java/snpEff/snpEff.jar eff -noStats $snpeff_genome -ud 50000 $input > $base.50kb.vcf

####START OF THE CHOOSE YOUR OWN ADVENTURE####
	echo
	while true
	do
	read -p "Want to go down the lectin or PGLYRP pathway? (type l or p) " answer
	case $answer in

#####LECTIN PATHWAY#####
	[lL]* )
			#variables for lectins or pglyrps
			family=lectin
			gene_interval=/media/russ/data/bovine/ref_files/lectin_intervals.bed
						
			echo "Lectins it is!"
			echo "Let's do it!"
			mkdir -p $family
			#MAKE VCF FILES WITH ONLY LECTINS
			java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R $gatk_genome -L $gene_interval -V $base.vcf -o $family/$base.all.$family.vcf
			java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R $gatk_genome -L $gene_interval -V $base.50kb.vcf -o $family/$base.all.$family.50kb.vcf
			
			#SPLIT INTO NORMAL AND DISEASED FOR BOTH UPSTREAM TYPES
			java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R $gatk_genome -sn group1 -sn group2 -sn group3 -sn group4 -sn group5 -sn group6 -sn group7 -sn group8 -V $family/$base.all.$family.vcf -o $family/$base.normal.$family.vcf
			java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R $gatk_genome -sn group9 -sn group10 -sn group11 -sn group12 -sn group13 -sn group14 -sn group15 -sn group16 -sn group17 -sn group18 -sn group19 -sn group20 -sn group21 -sn group22 -sn group23 -sn group24 -V $family/$base.all.$family.vcf -o $family/$base.diseased.$family.vcf
			java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R $gatk_genome -sn group1 -sn group2 -sn group3 -sn group4 -sn group5 -sn group6 -sn group7 -sn group8 -V $family/$base.all.$family.50kb.vcf -o $family/$base.normal.$family.50kb.vcf
			java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R $gatk_genome -sn group9 -sn group10 -sn group11 -sn group12 -sn group13 -sn group14 -sn group15 -sn group16 -sn group17 -sn group18 -sn group19 -sn group20 -sn group21 -sn group22 -sn group23 -sn group24 -V $family/$base.all.$family.50kb.vcf -o $family/$base.diseased.$family.50kb.vcf
			
			#MAKE SUBDIRECTORIES FOR ALL THE DIFFERENT REGIONS & GENES
			for x in coding intron downstream upstream 50kb_up 5utr 3utr
			do
				mkdir -p $family/$x
				for gene_array in colec10 colec11 colec12 colec43 colec46 cgn mbl1 mbl2 sftpa sftpd fcn1 masp1 masp2
				do
					mkdir -p $family/$x/$gene_array
				done
			done
			
			#TIDY UP
			mv $family/*50kb.vcf* $family/50kb_up
			
			#Files are now prepped and ready to be split into regions & genes.
			
			#Hash with ensembl references and gene names
			declare -A gene_hash
			gene_hash=(
				["ENSBTAG00000017343"]="colec10"
				["ENSBTAG00000016225"]="colec11"
				["ENSBTAG00000007705"]="colec12"
				["ENSBTAT00000066000.1"]="fcn1" #decided on this transcript variant with BL
				# FCN3 not found
				["ENSBTAG00000012467"]="masp1"
				["ENSBTAG00000012808"]="masp2"
				["ENSBTAT00000001165.3"]="mbl1" #transcript_reference - lots of confusion with SFTPA1 and MBl1
				["ENSBTAG00000007049"]="mbl2"
				["ENSBTAT00000031298.3"]="sftpa"	#transcript_reference
				["ENSBTAG00000046421"]="sftpd"
				["ENSBTAT00000018649.5"]="cgn" #two transcripts, decided on thsi one with BL
				["ENSBTAG00000048082"]="colec46"
				["ENSBTAG00000047317"]="colec43"
				)

			
			for nd in normal diseased all
			do
				for x in coding intron downstream upstream 50kb_up 5utr 3utr
				do
					for snpeff_gene in ${!gene_hash[@]}
					do
						gene_array="${gene_hash[$snpeff_gene]}"

						if [[ $gene_array == "mbl1" ]] || [[ $gene_array == "sftpa" ]] || [ $gene_array == "fcn1" ] || [ $gene_array == "cgn" ]
						then
							ann_type=ANN[*].TRID
						else
							ann_type=ANN[*].GENEID
						fi
# 					for gene_array in colec10 colec11 colec12 colec43 colec46 cgn mbl1 mbl2 sftpa sftpd fcn1 fcn3 masp1 masp2
# 					do
# 					
# 						ann_type=ANN[*].GENE
# 						if [[ $gene_array == "colec43" ]]
# 						then
# 							snpeff_gene=ENSBTAG00000047317
# 						elif [[ $gene_array == "colec46" ]]
# 						then
# 							snpeff_gene=ENSBTAG00000048082
# 						elif [[ $gene_array == "cgn" ]]
# 						then
# 							snpeff_gene=ENSBTAG00000006536
# 						elif [[ $gene_array == "mbl1" ]]
# 						then
# 							snpeff_gene=ENSBTAT00000001165.3
# 							ann_type=ANN[*].TRID
# 						elif [[ $gene_array == "sftpa" ]]
# 						then
# 							snpeff_gene=ENSBTAT00000031298.3
# 							ann_type=ANN[*].TRID
# 						elif [[ $gene_array == "fcn1" ]]
# 						then
# 							snpeff_gene=ENSBTAG00000048155
# 						elif [[ $gene_array == "fcn3" ]]
# 						then
# 							snpeff_gene=FIX_ME_PLEASE_THIS_ISNT_RIGHT
# 						else
# 							snpeff_gene=$(echo $gene_array | tr [:lower:] [:upper:])
# 						fi


						if [[ $x != "50kb_up" ]]
						then
							if [[ $x == "intron" ]]
							then
								type=intron_variant							
								echo $nd,$x,$gene_array, $type, $ann_type
								java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && ($ann_type = '$snpeff_gene')" $family/$base.$nd.$family.vcf > $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf
							elif [[ $x == "upstream" ]]
							then
								echo $nd,$x,$gene_array, $type, $ann_type
								type="$x"_gene_variant
								java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && ($ann_type = '$snpeff_gene') | (ANN[*].EFFECT has '5_prime_UTR_variant') && ($ann_type = '$snpeff_gene')" $family/$base.$nd.$family.vcf > $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf
							elif [[ $x == "5utr" ]]
							then
								echo $nd,$x,$gene_array, $type, $ann_type
								type="$x"_gene_variant
								java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '5_prime_UTR_variant') && ($ann_type = '$snpeff_gene')" $family/$base.$nd.$family.vcf > $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf
							elif [[ $x == "3utr" ]]
							then
								echo $nd,$x,$gene_array, $type, $ann_type
								type="$x"_gene_variant
								java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '3_prime_UTR_variant') && ($ann_type = '$snpeff_gene')" $family/$base.$nd.$family.vcf > $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf
							elif [[ $x == "downstream" ]]
							then
								echo $nd,$x,$gene_array, $type, $ann_type
								type="$x"_gene_variant
								java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && ($ann_type = '$snpeff_gene') | (ANN[*].EFFECT has '3_prime_UTR_variant') && ($ann_type = '$snpeff_gene')" $family/$base.$nd.$family.vcf > $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf
							else
								echo $nd,$x,$gene_array, $type, $ann_type
								java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') && ($ann_type = '$snpeff_gene') | (ANN[*].EFFECT has 'synonymous_variant') && ($ann_type = '$snpeff_gene') | (ANN[*].EFFECT has 'stop_gained') && ($ann_type = '$snpeff_gene')" $family/$base.$nd.$family.vcf > $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf
							fi
						else
							echo $nd,$x,$gene_array, $type, $ann_type
							java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'upstream_gene_variant') && ($ann_type = '$snpeff_gene') | (ANN[*].EFFECT has '5_prime_UTR_variant') && ($ann_type = '$snpeff_gene')" $family/50kb_up/$base.$nd.$family.50kb.vcf > $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf
						fi
					done
				done
			done

			break;;

######PGLYRP PATHWAY####
	[pP]* ) 
			echo "Ok, PGLYRPS. Kind of lame, but whatever"
			echo "PGLYRP-ing..."
			#variables for lectins or pglyrps
			family=pglyrp
			gene_interval=/media/russ/data/bovine/ref_files/pglyrp_intervals.intervals
						
			mkdir -p $family
			#MAKE VCF FILES WITH ONLY PGLYRPS OR LECTINS
			java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R $gatk_genome -L $gene_interval -V $base.vcf -o $family/$base.all.$family.vcf
			java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R $gatk_genome -L $gene_interval -V $base.50kb.vcf -o $family/$base.all.$family.50kb.vcf
			
			#SPLIT INTO NORMAL AND DISEASED
			java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R $gatk_genome -sn group1 -sn group2 -sn group3 -sn group4 -sn group5 -sn group6 -sn group7 -sn group8 -V $family/$base.all.$family.vcf -o $family/$base.normal.$family.vcf
			java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R $gatk_genome -sn group9 -sn group10 -sn group11 -sn group12 -sn group13 -sn group14 -sn group15 -sn group16 -sn group17 -sn group18 -sn group19 -sn group20 -sn group21 -sn group22 -sn group23 -sn group24 -V $family/$base.all.$family.vcf -o $family/$base.diseased.$family.vcf
			java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R $gatk_genome -sn group1 -sn group2 -sn group3 -sn group4 -sn group5 -sn group6 -sn group7 -sn group8 -V $family/$base.all.$family.50kb.vcf -o $family/$base.normal.$family.50kb.vcf
			java -jar ~/java/GATK/GenomeAnalysisTK.jar -T SelectVariants -R $gatk_genome -sn group9 -sn group10 -sn group11 -sn group12 -sn group13 -sn group14 -sn group15 -sn group16 -sn group17 -sn group18 -sn group19 -sn group20 -sn group21 -sn group22 -sn group23 -sn group24 -V $family/$base.all.$family.50kb.vcf -o $family/$base.diseased.$family.50kb.vcf
			
			for x in coding intron downstream upstream 50kb_up
			do
				mkdir -p $family/$x
				for gene_array in pglyrp1 pglyrp2 pglyrp3 pglyrp4
				do
					mkdir -p $family/$x/$gene_array
				done
			done
			
			mv $family/*50kb.vcf* $family/50kb_up
			
			#Files are now prepped and ready to be split into regions & genes.
# 			colec46=ENSBTAG00000048082
# 			cgn=ENSBTAG00000006536
# 			colec43=ENSBTAG00000047317
# 			mbl1=ENSBTAT00000001165 #transcript reference - gene reference is combined with SFTPA1 in ensembl
# 			sftpa=ENSBTAT00000031298 #transcript reference - gene reference is combined with SFTPA1 in ensembl
# 			fcn1=ENSBTAG00000048155
			
			for nd in normal diseased
			do
				for x in coding intron downstream upstream 50kb_up
				do
					for gene_array in pglyrp1 pglyrp2 pglyrp3 pglyrp4
					do
						if [[ $gene_array == "pglyrp4" ]]
						then
							snpeff_gene=
						else
							snpeff_gene=$(echo $gene_array | tr [:lower:] [:upper:])
						fi
						if [[ $x != "50kb_up" ]]
						then
							if [[ $x == "intron" ]]
							then
								type=intron_variant							
								echo sorting:$nd,$x,$gene_array
								java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$snpeff_gene')" $family/$base.$nd.$family.vcf > $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf
							elif [[ $x == "downstream" ]] || [[ $x == "upstream" ]]
							then
								echo sorting:$nd,$x,$gene_array
								type="$x"_gene_variant
								java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has '$type') && (ANN[*].GENE = '$snpeff_gene')" $family/$base.$nd.$family.vcf > $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf
							else
								echo sorting:$nd,$x,$gene_array
								java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') && (ANN[*].GENE = '$snpeff_gene') | (ANN[*].EFFECT has 'synonymous_variant') && (ANN[*].GENE = '$snpeff_gene') | (ANN[*].EFFECT has 'stop_gained') && (ANN[*].GENE = '$snpeff_gene')" $family/$base.$nd.$family.vcf > $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf
							fi
						else
							echo sorting:$nd,$x,$gene_array
							java -jar ~/java/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'upstream_gene_variant') && (ANN[*].GENE = '$snpeff_gene')" $family/50kb_up/$base.$nd.$family.50kb.vcf > $family/$x/$gene_array/$base.$nd.$x.$gene_array.vcf
						fi
					done
				done
			done			
			break;;
	* ) echo "You gotta enter l or p"
	esac
done

echo
echo
else
echo -n 
fi

read -p "Run vcf-compare (assumes snpeff has been run)? (y or Y) " vcf
	if [[ $vcf == "y" ]] || [[ $vcf == "yes" ]]
	then
		echo "Ok, running vcf-compare on everything."
		. vcf-compare_subscript.sh
	fi
	
read -p "How about summarizing all the snps? (y or Y) " summary
	if [[ $summary == "y" ]] || [[ $summary == "Y" ]]
	then
		echo "Great, summarizing."
		. snp_summarizer.sh
	else
		echo "k, bye"
	fi
fi

