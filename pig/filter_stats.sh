#Script to output the # of variants passing filter, the number of variants failing filters, and how many variants with rsIDs were found


#check if the summary table exists

#default directory
table=/media/russ/data/porcine/4_variants/groups/group_stats.txt

input=group.all.snps.vcf

#check if table already exists
if [ ! -f $table ]
then
	touch "$table"
	echo -e strategy'\t'total'\t'passed'\t'failed'\t'qd_filter'\t'mq_filter'\t'sort_filter'\t'fs_filter'\t'depth_filter'\t'dp-ac_filter'\t'readpos_filter'\t'with_rsid'\t'passed_rsid'\t'failed_rsid'\t'lowqual'\t'filtered_lowqual >> "$table"
fi

total=$(grep -cv \# "$input")
rsid=$(grep -v \# "$input" | grep -c rs)
lowqual=$(grep -v \# "$input" | grep -ic lowqual)


for file in q*vcf
do
	base="${file%.vcf}"

#gather stats
	echo working on "$base"...
# 	total=$(grep -v \# $file | wc -l)
	passed=$(grep -v \# $file | grep -ci pass)
	failed=$(grep -v \# $file | grep -vci pass)
	qdfilter=$(grep -v \# $file | grep -ci qdfilter)
	mqfilter=$(grep -v \# $file | grep -ci mqfilter)
	sortfilter=$(grep -v \# $file | grep -ci sorfilter)
	fsfilter=$(grep -v \# $file | grep -ci fsfilter)
	dpacfilter=$(grep -v \# $file | grep -ci dp-ac)
	readposfilter=$(grep -vc \# $file | grep -ci readposfilter)
	depth=$(grep -v \# $file | grep -ci depthfilter)

	#true positives and false negatives
	echo still going...
# 	rsid=$(grep -v \# $file | grep -c rs)
	rsidpassed=$(grep -v \# $file | grep rs | grep -ic pass)
	rsidfailed=$(grep -v \# $file | grep rs | grep -vci pass)
 
	#true negatives and false positives
	echo almost there...
# 	lowqual=$(grep -v \# $file | grep -ic lowqual)
	filteredlowqual=$(grep -v \# $file | grep -i lowqual | grep -ic filter)


	#pass stats to table
	echo -e "$base"'\t'"$total"'\t'"$passed"'\t'"$failed"'\t'"$qdfilter"'\t'"$mqfilter"'\t'"$sortfilter"'\t'"$fsfilter"'\t'"$depth"'\t'"$dpacfilter"'\t'"$readposfilter"'\t'"$rsid"'\t'"$rsidpassed"'\t'"$rsidfailed"'\t'"$lowqual"'\t'"$filteredlowqual" >> "$table"
	echo done!

done
