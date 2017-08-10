library(tidyverse)
library(viridis)
library(plyr)

annotation_info <- read.table("~/equine/2014_11_24/depth_of_regions/annotation_info.bed", h=T, sep="\t", na.strings = "na", stringsAsFactors = F)

#contains PGLYRPS
annotated_depth <- read.table("~/equine/2014_11_24/depth_of_regions/next_version.txt", h=T, sep="\t")

#no pglyrps and has fcn1-like.
annotated_depth <- read.table("~/equine/2014_11_24/depth_of_regions/annotated_depth_and_variants.txt", h=T, sep="\t")

#Clean up the annotation_info - there are a few records that have "na-" in them - didn't want to lose the info, but also don't want to use it. Check rows before proceeding... Don't actually have to do this step, but it will generate warnings down the line.
# 
# annotation_info$feature_start[annotation_info$gene_name == "SFTPD"][17:18] <- NA
# annotation_info$feature_end[annotation_info$gene_name == "SFTPD"][17:18] <- NA
# annotation_info$feature_start[annotation_info$gene_name == "FCN3"][17:18] <- NA
# annotation_info$feature_end[annotation_info$gene_name == "FCN3"][17:18] <- NA


#Function to summarize the data into a dataframe.
summarize <- function(info, data, depth){

    #Get a list of gene names and regions present within the annotation data frame
    genes <- unique(info$gene_name)
    regions <- unique(info$feature)

    #Create the output dataframe.
    output <- data.frame("Region" = as.character(), "Gene" = as.character(), "Total" = as.numeric(), "Sequenced" = as.numeric(), "No.variants" = as.numeric(), "Rate" = as.numeric(), stringsAsFactors = F)
    
    region_total <- 0
    
    #Ok I know you're not supposed to use for loops in R but this just seemed to accomplish what I needed in the time I needed to accomplish it... Plus no one will probably ever look at this, so great.
    for(region in regions){
      all_gene_targeted_total <- 0
      all_gene_sequenced_total <- 0
      all_gene_variants <- 0
      for(gene in genes){
        
        #Calculate the total targeted area for each gene and each region
        region_total <- sum((as.integer(subset(info, info$feature == region & info$gene_name == gene)$feature_end)- as.integer(subset(info, info$feature == region & info$gene_name == gene)$feature_start) + 1), na.rm=T)
        all_gene_targeted_total <- all_gene_targeted_total + region_total
        
        #Calculate the number of bases sequenced to a depth >= $depth (e.g. 5x / animal)
        region_sequenced <- nrow(subset(data, data$feature == region & data$gene_name == gene & data$depth >= depth))
        all_gene_sequenced_total <- all_gene_sequenced_total + region_sequenced
        
        #Calculate the number of variants for each gene/region
        variant_count <- sum(subset(data, data$feature == region & data$gene_name == gene)$is_variant)
        all_gene_variants <- all_gene_variants + variant_count
        
        #Calculate variation rate by gene/region
        rate <- variant_count/region_sequenced
        
        #Enter data into output dataframe
        output[nrow(output)+1,] <- c(region, gene, region_total, region_sequenced, variant_count, rate)
        
      }
      #Calculate the total targeted area by region for all genes together
      total_rate <- all_gene_variants/all_gene_sequenced_total
      output[nrow(output)+1,] <- c(region, "Total", all_gene_targeted_total, all_gene_sequenced_total, all_gene_variants, total_rate)
    }
  output$Total <- as.integer(output$Total)
  output$Sequenced <- as.integer(output$Sequenced)
  output$No.variants <- as.integer(output$No.variants)
  output$Rate <- as.numeric(output$Rate)
  return(output)
}

###I would argue that our variance rate will be elevated, because most estimates of the rate of variance are in the genome of one animal - our resolution maxes out at 3. So if you assume that some variants will be shared, and others not, then having multiple animals in a pool will result in a higher number of total variants, artifically elevating the rate. So, it is perhaps not useful to compare rates against other species, but rather between genes. So, for now, this is commented out and not pursued. Uncomment if you want to have a look at individual group stuff. 

#Can have a look at the different group variation rate by importing an file annotated by the variants present in each group (split the variants in GATK first, then run the add_variants_to_... python script using each sub-vcf file)

##Read in the group annotated files
# file.list <- list.files(path="~/equine/2014_11_24/depth_of_regions/", pattern="annotated.g")
# path="~/equine/2014_11_24/depth_of_regions/"
# all_files <- paste(path, file.list, sep="")
# group <- read.table("~/equine/2014_11_24/depth_of_regions/annotated.g5.txt", h=T, sep="\t")
##The data is read in as a list of dataframes, one for each condition. Labels -must- be applied using the python script "assign_genes_to_plink_output.py" before hand. The filenames are not used here to impart information to the data frames.
# all_groups <- lapply(all_files, read.table, h=T, sep="\t", stringsAsFactors=F)
# 
# b <- ldply(all_groups, summarize, info=annotation_info, depth=445)


##Apply the summarize function to the data
a <- summarize(annotation_info, annotated_depth, 445)
a

a_low <- summarize(annotation_info, annotated_depth, 10)
a_low

#Careful with subseting, the totals are NO LONGER ACCURATE!
colec <- subset(a, a$Gene != "PGLYRP1a" & a$Gene != "PGLYRP1b" & a$Gene != "PGLYRP1x" & a$Gene != "PGLYRP2" & a$Gene != "PGLYRP3" & a$Gene != "PGLYRP4")
pglyrp <- subset(a, a$Gene == "PGLYRP1a" | a$Gene == "PGLYRP1b" | a$Gene == "PGLYRP1x" | a$Gene == "PGLYRP2" | a$Gene == "PGLYRP3" | a$Gene == "PGLYRP4")




####STATS####
#Subset the data
#By Region:
exons <- subset(a, a$Region == "exon" & a$Gene != "Total")$Rate
introns <- subset(a, a$Region == "intron" & a$Gene != "Total")$Rate
downstream_3 <- subset(a, a$Region == "downstream_3" & a$Gene != "Total")$Rate
upstream_5 <- subset(a, a$Region == "upstream_5" & a$Gene != "Total")$Rate
upstream_50 <- subset(a, a$Region == "upstream_50" & a$Gene != "Total")$Rate

region_list <- list(exons, introns, downstream_3, upstream_5, upstream_50)

#By gene:
gene_list <- list()
for(gene in unique(a$Gene)){
  gene_list <- c(gene_list, list(subset(a, a$Gene == gene)$Rate))
}
gene_list

##Test each region for normality. low p-value indicates not normal.
lapply(region_list, shapiro.test)
lapply(gene_list, shapiro.test)

##Gene data is normal, but with only 5 observations is that reliable? 
aov_genes <- aov(a$Rate ~ a$Gene)
summary(aov_genes)
TukeyHSD(aov_genes)

##Region data is not normal. Test between regions can't use ANOVA, use non-parametric Kruskal-Wallis alternative.
kruskal.test(region_list)



##Output tables
write.table(b, file="~/Desktop/table.txt", sep="\t", row.names = F, quote = F)


##Plots
ggplot(subset(a, a$Gene != "Total"), aes(x=Gene)) + geom_bar(aes(y=Total, alpha = 0.6, fill = Region), stat="identity") + geom_bar(aes(y=Sequenced, fill = Region), stat="identity") + theme(axis.text.x = element_text(angle=90))
str(a)

ggplot(b, aes(x=Region)) + geom_boxplot(aes(y=Rate)) + ylim(c(0,0.05))
ggplot(b, aes(x=Region, y=Rate, group=Gene)) + geom_bar(aes(fill = Gene), stat="identity", position="dodge") + theme(axis.text.x = element_text(angle=90)) + scale_colour_viridis(discrete = T) + theme_bw()

##Collectins only
ggplot(colec, aes(x = Gene, y = Rate)) + 
  geom_bar(aes(fill = Region), stat="identity") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggplot(colec, aes(x = Region, y = Rate)) + 
  geom_boxplot() +
  theme_bw() + 
  ylab("Rate (SNP/bp)") + 
  xlab("")
