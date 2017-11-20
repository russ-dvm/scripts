library(tidyverse)
library(plyr)

# Read in all the snpEff results

file.list <- list.files(path="~/bovine/merged_runs/6_variants/lectin/byGene/", pattern="stats.genes.txt", recursive = T)
path="~/bovine/merged_runs/6_variants/lectin/byGene/"
all_files <- paste(path, file.list, sep="")

## Iterate over the list of files and read them all into a DF. If results are a bit wonky remember that the header in the original files contains a "#" - can remove using something like for d in */; do echo $d; base="${d%/}"; sed -i 's/\#G/G/' "$d""$base"*genes*; done
allData <- lapply(all_files, read.table, h=T)

## Combine each file into a single data frame. Will introduce NAs for categories that are present in one gene but not another.
main <- ldply(allData, bind_rows)



## Annotate missing gene names and remove extraneous genes (from 50kb region)
main$GeneName <- gsub("ENSBTAG00000006536", "CGN1", main$GeneName)
main$GeneName <- gsub("ENSBTAG00000047317", "COLEC43", main$GeneName)
main$GeneName <- gsub("ENSBTAG00000048082", "COLEC46", main$GeneName)
main$GeneName <- gsub("ENSBTAG00000048155", "FCN1", main$GeneName)
main$GeneName <- gsub("ENSBTAG00000046421", "SFTPD", main$GeneName)
main[grep("ENSBTAT00000001165", main$TranscriptId),]$GeneName <- "MBL1"

main <- main[-grep("RPS7", main$GeneName),]
main <- main[-grep("CD164L2", main$GeneName),]
main <- main[-grep("TARDBP", main$GeneName),]
main <- main[-grep("SRM", main$GeneName),]
main <- main[-grep("MAP3K6", main$GeneName),]




## Remove impact descriptors
impact <- grep("impact", colnames(main))
main <- main[, -impact]

## Create a list of the column names
a<-colnames(main)[5:ncol(main)]

## Tidy the data
main_gathered <- gather(main, type, count, a)
main_gathered[is.na(main_gathered$count),]$count <- 0

## Clean the category names
main_gathered$type <- gsub("variants_effect_", "", main_gathered$type)
main_gathered$type <- gsub("_variant", "", main_gathered$type)
main_gathered$type <- gsub("_", " ", main_gathered$type)

## Plot
ggplot(main_gathered, aes(x = GeneName, y = count)) + geom_bar(aes(fill = type), stat="identity") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1))


