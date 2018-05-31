library(tidyverse)
library(plyr)
library(gridExtra)
library(ggpubr)

all <- read.table("~/bovine/merged_runs/depth/depth.minqual.annotated.variants.txt", h = T, sep = "\t", stringsAsFactors = F)
annotation_info <- read.table("~/bovine/merged_runs/depth/cow_annotation.bed", h = T, sep = "\t", stringsAsFactors = F)


#Function to summarize the data into a dataframe.
summarize <- function(info, var_data, depth_cutoff){
  
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
      region_total <- sum((as.integer(subset(info, info$feature == region & info$gene_name == gene)$end)- as.integer(subset(info, info$feature == region & info$gene_name == gene)$start) + 1), na.rm=T)
      all_gene_targeted_total <- all_gene_targeted_total + region_total
      
      #Calculate the number of bases sequenced to a depth >= $depth (e.g. 5x / animal)
      region_sequenced <- nrow(subset(var_data, var_data$feature == region & var_data$gene_name == gene & var_data$depth >= depth_cutoff))
      all_gene_sequenced_total <- all_gene_sequenced_total + region_sequenced
      
      #Calculate the number of variants for each gene/region
      variant_count <- sum(subset(var_data, var_data$feature == region & var_data$gene_name == gene)$is_variant)
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

colec <- summarize(annotation_info, all, 600)
colec$rateKb <- colec$Rate*1000



#### STATS ####

## Create lists of the data by region and by gene
## Some excess data in the OG dataframe - trim out redundancies (eg 5UTR, 3UTR, 3kb_from_gene)

gene_list <- list()
for(gene in unique(colec$Gene)) {
  if (gene == "Total"){}
  else {
    colecTrimmed <- subset(colec, Region != "3kb_from_gene")
    colecTrimmed <- subset(colecTrimmed, Region != "3UTR")
    colecTrimmed <- subset(colecTrimmed, Region != "5UTR")
    colecTrimmed <- subset(colecTrimmed, Region != "5kb_from_gene")
    colecTrimmed <- subset(colecTrimmed, Region != "45kb_from_gene")
    tmp <- list(subset(colecTrimmed, colecTrimmed$Gene == gene)$Rate)
    names(tmp) <- gene
    gene_list <- c(gene_list, tmp)
  }
}

colecTrimmedByGene <- colecTrimmed
colecTrimmed <- subset(colecTrimmed, Gene != "Total")

region_list <- list()
for(region in unique(colecTrimmed$Region)) {
  tmp <- list(subset(colecTrimmed, colecTrimmed$Region == region)$Rate)
  names(tmp) <- region
  region_list <- c(region_list, tmp)
}

lapply(region_list, shapiro.test)
lapply(gene_list, shapiro.test)

#Use kruskal wallis - low sample numbers and questionable normality
kruskal.test(region_list)
kruskal.test(gene_list)

#Dunn test for post-hoc evaluation
library(dunn.test)

#Irritatingly this doens't use the list names in the output...
dunn.test(gene_list, label = T, list = T, method = "bh")

#Make a dataframe and gather it to solve the issue
geneDunn <- ldply(gene_list) %>% gather(key = ".id")
colnames(geneDunn) <- c("gene", "region", "value")
dunnRes <- dunn.test(x = geneDunn$value, g = geneDunn$gene, method = "bh", list = T)
sigComps <- dunnRes$comparisons[dunnRes$P.adjusted<0.05]
sigCompsforPubr <- strsplit(sigComps, " - ")

#Manual reordering of the significant comparisons for a more aesthetically pleasing plot... will have to be altered if the data changes.
sigRe <-c(sigCompsforPubr[3], sigCompsforPubr[4], sigCompsforPubr[8], sigCompsforPubr[5], sigCompsforPubr[6], sigCompsforPubr[7], sigCompsforPubr[1], sigCompsforPubr[2])

## DFs for plotting - get rid of totals - one with UTRs and one wihtout
colecUTR <- subset(colec, Gene != "Total")
colecUTR <- subset(colecUTR, Region != "5kb_from_gene")
colecUTR <- subset(colecUTR, Region != "45kb_from_gene")
colecUTR <- subset(colecUTR, Region != "3kb_from_gene")
colecUTR$Region <- gsub("5UTR", "5' UTR", colecUTR$Region)
colecUTR$Region <- gsub("3UTR", "3' UTR", colecUTR$Region)

## Relevel the regions
colecTrimmed$Region <- gsub("45kb_from_start", "5-50 kb upstream", colecTrimmed$Region)
colecTrimmed$Region <- gsub("5kb_from_start", "Upstream 5 kb", colecTrimmed$Region)
colecTrimmed$Region <- gsub("exon", "Coding", colecTrimmed$Region)
colecTrimmed$Region <- gsub("intron", "Intron", colecTrimmed$Region)
colecTrimmed$Region <- gsub("3kb_from_stop", "Downstream 3 kb", colecTrimmed$Region)
colecTrimmed$Region <- factor(colecTrimmed$Region, levels = c("5-50 kb upstream", "Upstream 5 kb", "Coding", "Intron", "Downstream 3 kb"))

####PLOTS####
a <- ggplot(colecTrimmed, aes(x = reorder(colecTrimmed$Gene, colecTrimmed$Rate, function(x)-median(x)), y = Rate*1000)) + 
  geom_boxplot() + 
  # stat_summary(fun.y = mean, colour = "grey", geom = "text", label = "----") +
  theme_classic() +
  theme(text = element_text(size = 10)) +
    ylab("Variant density (per kb) by gene") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust = -0.14, size = 10)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) +
  annotate("text", label = paste(expression("\u2020")), x = 1, y = 25) +
  annotate("text", label = paste(expression("\u2021")), x = 2, y = 25) + 
  annotate("text", label = paste(expression("\u2020")), x = 3, y = 25) +
  stat_summary(fun.y = mean, colour = "grey56", geom = "text", label = "....", size = 4)
  
a

b<- ggplot(subset(colecTrimmed, colecTrimmed$Gene != "Total"), aes(x = Region, y = Rate*1000)) +
  geom_boxplot() +
  theme_classic() +
  ylab("Variant density (per kb) by region") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle("b)") + 
  theme(plot.title = element_text(hjust = -0.14, size = 10)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) +
  stat_summary(fun.y = mean, colour = "grey56", geom = "text", label = "..........", size = 4)

b

grid.arrange(a, b, ncol = 2)

## Convoluted way to get the plots to be the same height (gridarrange was not doing the trick) ** ONLY needed if both boxplots are used. the Bar plot is no problemo.
## convert plots to gtable objects
library(gtable)
library(grid) # low-level grid functions are required
g1 <- ggplotGrob(a)
g2 <- ggplotGrob(b)
g <- cbind(g1, g2, size="first") # stack the two plots
# g$heights <- unit.pmin(g1$heights, g2$heights) # use the smallest heights

fig_size <- c(190, 120) # inches
margin <- unit(1, "line")

g$vp <- viewport(width = unit(fig_size[1], "mm"), height=unit(fig_size[2],"mm")- margin)

library("Cairo")
ggsave("~/Dropbox/temp/figure_RATES.pdf", g, width=fig_size[1], height=fig_size[2], units = "mm", dpi = 1000, device = cairo_pdf)


### Write out the table
## Relevel the regions
colecTrimmedByGene$Region <- gsub("45kb_from_start", "5-50 kb upstream", colecTrimmedByGene$Region)
colecTrimmedByGene$Region <- gsub("5kb_from_start", "Upstream 5 kb", colecTrimmedByGene$Region)
colecTrimmedByGene$Region <- gsub("exon", "Coding", colecTrimmedByGene$Region)
colecTrimmedByGene$Region <- gsub("intron", "Intron", colecTrimmedByGene$Region)
colecTrimmedByGene$Region <- gsub("3kb_from_stop", "Downstream 3 kb", colecTrimmedByGene$Region)
colecTrimmedByGene$Region <- factor(colecTrimmedByGene$Region, levels = c("5-50 kb upstream", "Upstream 5 kb", "Coding", "Intron", "Downstream 3 kb"))
colecTrimmedByGene <- colecTrimmedByGene[order(colecTrimmedByGene$Gene),]
write.table(colecTrimmedByGene, file = "~/Dropbox/temp/rates.txt", row.names = F, quote = F, sep = "\t")

######SANKEY PLOT######
###Note: this doesn't work on the Unix comp - the github version of ggforce wouldn't install and I can't be bothered to resolve the dependencies - but does work on windows. Things to improve - vertically center the axes, adjust labels so that they appear outside of boxes? Too small to read ATM. Easy enough to add a third column for allele frequency significance.git c
# test <- gather_set_data(colecTrimmed, 1:2)
# ggplot(test, aes(x = x, id = id, value = No.variants, split = y)) + geom_parallel_sets(aes(fill = test$Gene), alpha = 0.3, show.legend = F) + geom_parallel_sets_axes(axis.width = 0.1) + geom_parallel_sets_labels(colour = "white", size = 2, angle = 0) + theme_bw() + theme(panel.grid = element_blank()) + xlab("") + ylab("") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())