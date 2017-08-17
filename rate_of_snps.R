library(tidyverse)
library(viridis)
library(plyr)
library(data.table)

annotation_info <- read.table("~/equine/2014_11_24/depth_of_regions/annotation_info.bed", h=T, sep="\t", na.strings = "na", stringsAsFactors = F)
annotation_info <- read.table("~/Desktop/annotation_info_utr_adjusted.bed", h=T, sep="\t", na.strings = "na", stringsAsFactors = F, quote = "")

#contains PGLYRPS
annotated_depth <- read.table("~/equine/2014_11_24/depth_of_regions/next_version.txt", h=T, sep="\t")

#no pglyrps and has fcn1-like.
annotated_depth <- read.table("~/Desktop/results_utr_adjusted.txt", h=T, sep="\t")

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

###I would argue that our variance rate will be elevated, because most estimates of the rate of variance are in the genome of one animal - our resolution maxes out at 3. So if you assume that some variants will be shared, and others not, then having multiple animals in a pool will result in a higher number of total variants, artifically elevating the rate. So, it is perhaps not useful to compare rates against other species, but rather between genes. 


##Apply the summarize function to the data
a <- summarize(annotation_info, annotated_depth, 445)
a

## Adjust for MBL1 50kb, which overlaps with SFTPD and it's upstream region. At some point when time is abundant, maybe integrate this into the function... 
mbl_50_start <- 88892827
mbl_50_end <- 88928053
a <- data.table(a)
a[Region == "upstream_50" & Gene == "MBL1"]$Total <- mbl_50_end - mbl_50_start
mbl50_sequenced <- nrow(subset(annotated_depth, pos >= mbl_50_start & pos <= mbl_50_end))
a[Region == "upstream_50" & Gene == "MBL1"]$Sequenced <- mbl50_sequenced
mbl50_var <- nrow(subset(annotated_depth, pos >= mbl_50_start & pos <= mbl_50_end & is_variant == T))
a[Region == "upstream_50" & Gene == "MBL1"]$No.variants <- mbl50_var
a[Region == "upstream_50" & Gene == "MBL1"]$Rate <- a[Region == "upstream_50" & Gene == "MBL1"]$No.variants/a[Region == "upstream_50" & Gene == "MBL1"]$Sequenced
##Adjust the upstream_50 total
a[Region == "upstream_50" & Gene == "Total"]$No.variants <- sum(a[Region == "upstream_50" & Gene != "Total"]$No.variants)

## The damn 5kb upstream is also off, since priority in the annotation file was given to the downstream 3 kb of SFTPA
annotation_info <- data.table(annotation_info)
mbl_5_end <- as.integer(annotation_info[feature == "upstream_5" & gene_name == "MBL1",]$feature_start) -1
##For the start pos, take the "end" of the original upstream_5 and substract 5000 to get what it should've been
mbl_5_start <- as.integer(annotation_info[feature == "upstream_5" & gene_name == "MBL1",]$feature_end) - 5000
mbl_5_adjust <- mbl_5_end - mbl_5_start

a[Region == "upstream_5" & Gene == "MBL1"]$Total <- a[Region == "upstream_5" & Gene == "MBL1"]$Total + mbl_5_adjust + 1
mbl5_sequenced <- nrow(subset(annotated_depth, pos >= mbl_5_start & pos <= mbl_5_end))
a[Region == "upstream_5" & Gene == "MBL1"]$Sequenced <- a[Region == "upstream_5" & Gene == "MBL1"]$Sequenced + mbl5_sequenced
mbl5_var <- nrow(subset(annotated_depth, pos >= mbl_5_start & pos <= mbl_5_end & is_variant == T))
a[Region == "upstream_5" & Gene == "MBL1"]$No.variants <- a[Region == "upstream_5" & Gene == "MBL1"]$No.variants + mbl5_var
a[Region == "upstream_5" & Gene == "MBL1"]$Rate <- a[Region == "upstream_5" & Gene == "MBL1"]$No.variants/a[Region == "upstream_5" & Gene == "MBL1"]$Sequenced


###Careful with subseting, the totals are NO LONGER ACCURATE!
colec <- subset(a, a$Gene != "PGLYRP1a" & a$Gene != "PGLYRP1b" & a$Gene != "PGLYRP1x" & a$Gene != "PGLYRP2" & a$Gene != "PGLYRP3" & a$Gene != "PGLYRP4")
# pglyrp <- subset(a, a$Gene == "PGLYRP1a" | a$Gene == "PGLYRP1b" | a$Gene == "PGLYRP1x" | a$Gene == "PGLYRP2" | a$Gene == "PGLYRP3" | a$Gene == "PGLYRP4")


##Revalue the factors to something more interpretable (requires plyr)
colec$Region <- revalue(colec$Region, c("downstream_3" = "Downstream 3 kb", "exon" = "Coding", "intron" = "Introns", "upstream_5" = "Upstream 5 kb", "upstream_50" = "Upstream 5-50 kb"))
colec$Region <- factor(colec$Region, levels = c("Upstream 5-50 kb", "Upstream 5 kb", "Coding", "Introns", "Downstream 3 kb"))

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


##Test each region for normality. low p-value indicates not normal.
lapply(region_list, shapiro.test)
lapply(gene_list, shapiro.test)

##Gene data is normal, but with only 5 observations is that reliable? 
aov_genes <- aov(a$Rate ~ a$Gene)
summary(aov_genes)
TukeyHSD(aov_genes)

##Region data is not normal. Test between regions can't use ANOVA, use non-parametric Kruskal-Wallis alternative.
kruskal.test(region_list)


##Plots
## Bar plot of overall variation by gene
ggplot(subset(a, a$Gene != "Total"), aes(x=Gene)) + geom_bar(aes(y=Total, alpha = 0.6, fill = Region), stat="identity") + geom_bar(aes(y=Sequenced, fill = Region), stat="identity") + theme(axis.text.x = element_text(angle=90))


##Collectins only - box and bar plots, by region/gene
ggplot(subset(colec, colec$Gene != "Total"), aes(x = Gene, y = Rate*1000)) + 
  geom_bar(aes(fill = Region), stat="identity", position = "dodge") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylab("Rate (SNPs/kb)")

ggplot(subset(colec, colec$Gene != "Total"), aes(x = Region, y = Rate)) + 
  geom_boxplot() +
  geom_jitter(aes(colour = Gene)) +
  theme_bw() +
  ylab("Rate (SNP/bp)") + 
  xlab("") + 
  scale_colour_viridis(discrete = T) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggplot(subset(colec, colec$Gene != "Total"), aes(x = Gene, y = Rate)) + 
  geom_boxplot() + 
  geom_jitter(aes(colour = Region)) +
  geom_line(aes()) +
  theme_bw() +
  ylab("Rate (SNP/bp)") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

  