library(tidyverse)
library(plyr)
library(data.table)
library(ggpubr)
library(gridExtra)

annotation_info <- read.table("~/equine/2014_11_24/depth_of_regions/annotation_info_utr_adjusted.bed", h=T, sep="\t", na.strings = "na", stringsAsFactors = F, quote = "")
# annotation_info <- read.table("~/Desktop/annotation_info_utr_adjusted.bed", h=T, sep="\t", na.strings = "na", stringsAsFactors = F, quote = "")

#contains PGLYRPS
# annotated_depth <- read.table("~/equine/2014_11_24/depth_of_regions/next_version.txt", h=T, sep="\t")

#no pglyrps and has fcn1-like.
annotated_depth <- read.table("~/equine/2014_11_24/depth_of_regions/results_utr_adjusted.txt", h=T, sep="\t")
annotated_depth <- read.table("~/equine/2014_11_24/depth_of_regions/results_utr_adjusted_aug30.txt", h=T, sep="\t")

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


####ADJUSTMENTS####
### Adjust for MBL1 50kb, which overlaps with SFTPD and it's upstream region. At some point when time is abundant, maybe integrate this into the function... 
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
colec <- subset(a, a$Gene != "PGLYRP1a" & a$Gene != "PGLYRP1b" & a$Gene != "PGLYRP1x" & a$Gene != "PGLYRP2" & a$Gene != "PGLYRP3" & a$Gene != "PGLYRP4" & a$Gene != "Total")
# pglyrp <- subset(a, a$Gene == "PGLYRP1a" | a$Gene == "PGLYRP1b" | a$Gene == "PGLYRP1x" | a$Gene == "PGLYRP2" | a$Gene == "PGLYRP3" | a$Gene == "PGLYRP4")




####STATS####
#Subset the data
#By Region:

exons <- subset(colec, colec$Region == "exon" & colec$Gene != "Total")$Rate
introns <- subset(colec, colec$Region == "intron" & colec$Gene != "Total")$Rate
downstream_3 <- subset(colec, colec$Region == "downstream_3" & colec$Gene != "Total")$Rate
upstream_5 <- subset(colec, colec$Region == "upstream_5" & colec$Gene != "Total")$Rate
upstream_50 <- subset(colec, colec$Region == "upstream_50" & colec$Gene != "Total")$Rate

region_list <- list("exons" = exons, "introns" = introns, "downstream_3" = downstream_3, "upstream_5" = upstream_5, "upstream_50" = upstream_50)

#By gene:
gene_list <- list()
for(gene in unique(colec$Gene)){
  if (gene == "Total"){}
  else{
  gene_list <- c(gene_list, list(subset(colec, Gene == gene)$Rate))
}}
names(gene_list) <- unique(colec$Gene)[c(1:12)]

gene_list[["FCN3"]][5] <- 0
gene_avg <- lapply(gene_list, mean)
gene_avg_df <- ldply(gene_avg, data.frame)
colnames(gene_avg_df) <- c("gene", "avg_var")
gene_avg_df$avg_var_per_kb <- gene_avg_df$avg_var * 1000

##Test each region for normality. low p-value indicates not normal.
lapply(region_list, shapiro.test)
lapply(gene_list, shapiro.test)

##Gene data is normal, but with only 5 observations is that reliable? 
aov_genes <- aov(a$Rate ~ a$Gene)
summary(aov_genes)
tuk_genes <- TukeyHSD(aov_genes)
plot(TukeyHSD(aov_genes))

tuk_table <- data.table(tuk_genes[[1]], keep.rownames = T)
colnames(tuk_table) <- c("rn", "diff", "lwr", "upr", "p")
tuk_table[p<0.05]

##Region data is not normal. Test between regions shouldn't use ANOVA, use non-parametric Kruskal-Wallis alternative.
region_list[[5]][2] <- 0

kruskal.test(region_list)

## AOV for brandon 
aov_region <- aov(a$Rate ~ a$Region)
summary(aov_region)
tuk_regions <- TukeyHSD(aov_region)
tuk_regions
plot(tuk_regions)

lapply(region_list, mean)


####PLOTS####

##Revalue the factors to something more interpretable (requires plyr). Makes the graphs nicer
colec$Region <- revalue(colec$Region, c("downstream_3" = "Downstream 3 kb", "exon" = "Coding", "intron" = "Introns", "upstream_5" = "Upstream 5 kb", "upstream_50" = "Upstream 5-50 kb"))
colec$Region <- factor(colec$Region, levels = c("Upstream 5-50 kb", "Upstream 5 kb", "Coding", "Introns", "Downstream 3 kb"))

##Fix the NA value 
colec[Gene == "FCN3" & Sequenced == 0]$Rate <- 0

## Bar plot of overall variation by gene
# ggplot(subset(a, a$Gene != "Total"), aes(x=Gene)) + geom_bar(aes(y=Total, alpha = 0.6, fill = Region), stat="identity") + geom_bar(aes(y=Sequenced, fill = Region), stat="identity") + theme(axis.text.x = element_text(angle=90))


## Collectins only - box and bar plots, by region/gene
# 
# ggplot(subset(colec, colec$Gene != "Total"), aes(x = Region, y = Rate*1000)) +
#   geom_bar(aes(fill = Gene), stat="identity", position = "dodge") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#   ylab("Number of variants per kb") +
#   scale_fill_grey() +
#   # theme(legend.position = c(1,1), legend.justification= c(1,1)) +
#   theme(legend.background = element_rect(fill = "white", linetype = "solid", colour = "black", size = 0.3)) +
#   theme(legend.text = element_text(face = "italic", size = 7), legend.title = element_text(size = 7)) +
#   theme(legend.position = c(1,1), legend.justification= c(1,1))
#   
# ggplot(subset(colec, colec$Gene != "Total"), aes(x = Gene, y = Rate*1000)) + 
#   geom_bar(aes(fill = Region), stat="identity", position = "dodge") +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#   ylab("Number of variants per kb") +
#   xlab("") +
#   scale_fill_grey() +
#   theme(legend.position = c(1,1), legend.justification= c(1,1)) +
#   theme(legend.text = element_text(size = 7), legend.title = element_text(size = 7)) + 
#   theme(legend.background = element_rect(fill = "white", linetype = "solid", colour = "black", size = 0.3)) +
#   ggtitle("b)") + 
#   theme(plot.title = element_text(hjust = -0.14, size = 10))



a <- ggplot(subset(colec, colec$Gene != "Total"), aes(x = Gene, y = Rate*1000)) + 
  geom_boxplot() + 
  theme(text = element_text(size = 10)) +
  stat_summary(fun.y = mean, colour = "grey", geom = "text", label = "----") +
  theme_bw() +
  ylab("Number of variants per kb") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  annotate("text", label = "*", x = 4, y = 32.5, size = 6) +
  annotate("text", label = "*", x = 5, y = 36, size = 6) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust = -0.14, size = 10)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank())
a

b<- ggplot(subset(colec, colec$Gene != "Total"), aes(x = Region, y = Rate*1000)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Number of variants per kb") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  stat_summary(fun.y = mean, colour = "grey", geom = "text", label = "----", size = 8) +
  ggtitle("b)") + 
  theme(plot.title = element_text(hjust = -0.14, size = 10)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank())


grid.arrange(a, b, ncol = 2)

## Convoluted way to get the plots to be the same height (gridarrange was not doin the trick) ** ONLY needed if both boxplots are used. the Bar plot is no problemo.
## convert plots to gtable objects
library(gtable)
library(grid) # low-level grid functions are required
g1 <- ggplotGrob(a)
g2 <- ggplotGrob(b)
g <- cbind(g1, g2, size="first") # stack the two plots
# g$heights <- unit.pmin(g1$heights, g2$heights) # use the smallest heights

fig_size <- c(190, 100) # inches
margin <- unit(1, "line")

g$vp <- viewport(width = unit(fig_size[1], "mm"), height=unit(fig_size[2],"mm")- margin)

library("Cairo")
ggsave("~/Dropbox/temp/figure_2.pdf", g, width=fig_size[1], height=fig_size[2], units = "mm", dpi = 1000, device = cairo_pdf)


####CHECKING GC CONTENT####
gc <- read.table("~/equine/2014_11_24/gc_content/gene_gc/summary.gc", sep = "\t")
colnames(gc) <- c("gene", "pct_gc")
tstv <- read.table("~/equine/2014_11_24/gc_content/gene_tstv/summary.tstv")
colnames(tstv) <- c("gene", "tstv")

merged <- merge(gc, tstv)
merged_sorted <- merged[order(merged$pct_gc),]

ggplot(merged_sorted, aes(x=tstv, y = pct_gc)) + geom_point() + geom_smooth(method = lm)
cor(merged$gc, merged$tstv)

sup5 <- gc
colnames(sup5) <- c("Target Region", "Percent GC")
sup5$`Target Region` <- c("COLEC10", "COLEC11", "COLEC12", "SFTPA, MBL1, SFTPD", "FCN3", "FCN1, FCN1-LIKE", "MASP1", "MASP2", "MBL2")

####INDEL VS SNP####
ind_snp <- read.table("~/equine/2014_11_24/gc_content/gene_vcfs/summary.txt", sep = "\t")
ind_snp_spread <- spread(ind_snp, V2, V3)
ind_snp_spread$ratio <- ind_snp_spread$indel/ind_snp_spread$snp
ggplot(ind_snp_spread, aes(x = indel, y = snp)) + geom_point(aes(colour = V1, shape = V1, size = 4)) + geom_smooth(method = lm) +
  scale_shape_manual(values = c(3:19))

cor(ind_snp_spread$snp, ind_snp_spread$indel)
## Strong correlation between SNPs and INDELs, which supports the thoery that there are more point mutations in areas that have more INDELs... 


## SUPP TABLE 5
sup5
colnames(ind_snp_spread) <- c("Target Region", "INDEL", "SNV", "Ratio")
ind_snp_spread$`Target Region` <-  c("COLEC10", "COLEC11", "COLEC12", "SFTPA, MBL1, SFTPD", "FCN1, FCN1-LIKE", "FCN3", "MASP1", "MASP2", "MBL2")
ind_snp_spread
sup5 <- merge(sup5, ind_snp_spread)
sup5
write.table(sup5, "~/Dropbox/temp/sup5.txt", row.names=F, quote=F, sep="\t")
