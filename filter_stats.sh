#Script to output the # of variants passing filter, the number of variants failing filters, and how many variants with rsIDs were found


#check if the summary table exists

#default directory
##PIG
# table=/media/russ/data/porcine/4_variants/filter_stats.txt
##COW
table=/media/russ/data/bovine/merged_runs/5_filter_testing/filter_stats.txt

#check if table already exists
if [ ! -f $table ]
then
	touch "$table"
	echo -e strategy'\t'total'\t'passed'\t'failed'\t'qd_filter'\t'mq_filter'\t'sort_filter'\t'fs_filter'\t'dp-ac_filter'\t'readpos_filter'\t'with_rsid'\t'passed_rsid'\t'failed_rsid'\t'lowqual'\t'filtered_lowqual >> "$table"
fi

#gather stats
echo working on it...
total=$(grep -v \# $1 | wc -l)
passed=$(grep -v \# $1 | grep -vci filter)
failed=$(grep -v \# $1 | grep -ci filter)
qdfilter=$(grep -v \# $1 | grep -ci qdfilter)
mqfilter=$(grep -v \# $1 | grep -ci mqfilter)
sortfilter=$(grep -v \# $1 | grep -ci sorfilter)
fsfilter=$(grep -v \# $1 | grep -ci fsfilter)
dpacfilter=$(grep -v \# $1 | grep -ci dp-ac)
readposfilter=$(grep -vc \# $1 | grep -ci readposfilter)

#true positives and false negatives
echo still going...
rsid=$(grep -v \# $1 | grep -c rs)
rsidpassed=$(grep -v \# $1 | grep rs | grep -ic pass)
rsidfailed=$(grep -v \# $1 | grep rs | grep -vci pass)
 
#true negatives and false positives
echo almost there...
lowqual=$(grep -v \# $1 | grep -ic lowqual)
filteredlowqual=$(grep -v \# $1 | grep -i lowqual | grep -ic filter)

#pass stats to table
echo -e "$1"'\t'"$total"'\t'"$passed"'\t'"$failed"'\t'"$qdfilter"'\t'"$mqfilter"'\t'"$sortfilter"'\t'"$fsfilter"'\t'"$dpacfilter"'\t'"$readposfilter"'\t'"$rsid"'\t'"$rsidpassed"'\t'"$rsidfailed"'\t'"$lowqual"'\t'"$filteredlowqual" >> "$table"
echo done!