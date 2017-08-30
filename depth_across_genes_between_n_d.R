library(tidyverse)

#######FILE IMPORT
depth <- read.table("~/equine/2014_11_24/depth_of_regions/all.txt.sample_interval_summary", h=T, sep="\t", stringsAsFactors = F)

#######CLEAN DATA
#Take only the mean counts
mean_depths <- depth[grep("mean", colnames(depth))]
mean_depths$target <- depth$Target
colnames(mean_depths) <- gsub("_mean_cvg", "", colnames(mean_depths))

##Add gene names
mean_depths$gene <- NA
mean_depths[grep("chr15:88547816-88602733", mean_depths$target),]$gene <- "COLEC11"
mean_depths[grep("chr19:25141928-25256816", mean_depths$target),]$gene <- "MASP1"
mean_depths[grep("chr1:43242494-43299035", mean_depths$target),]$gene <- "MBL2"
mean_depths[grep("chr1:88892827-89006558", mean_depths$target),]$gene <- "SFTPA_MBL1_SFTPD"
mean_depths[grep("chr25:36785097-36875451", mean_depths$target),]$gene <- "FCN1_FCN1-like"
mean_depths[grep("chr2:28409278-28419936", mean_depths$target),]$gene <- "FCN3"
mean_depths[grep("chr2:40412424-40446898", mean_depths$target),]$gene <- "MASP2"
mean_depths[grep("chr8:40892387-40987322", mean_depths$target),]$gene <- "COLEC12"
mean_depths[grep("chr9:62213522-62301779", mean_depths$target),]$gene <- "COLEC10"

mean_depths_gathered <- gather(mean_depths, group, avg, -target, -gene)


#######ASSIGN DISEASE STATUS
##Manually curated vector of diseased groups
diseased <- c("group7", "group8", "group9", "group10", "group11", "group12", "group13", "group20", "group18")

for (i in 1:nrow(mean_depths_gathered)){
  gname <- mean_depths_gathered$group[i]
  if (gname %in% diseased){
    mean_depths_gathered$status[i] <- "diseased"
  }
  else{
    mean_depths_gathered$status[i] <- "healthy"
  }
}

#####STATS
##Differences between genes
depth_aov <- aov(mean_depths_gathered$avg ~ mean_depths_gathered$gene)
summary(depth_aov)
TukeyHSD(depth_aov)

##Differences between status at each gene
healthy <- subset(mean_depths_gathered, status == "healthy")
healthy_list <- split(healthy$avg, healthy$gene)
dis <- subset(mean_depths_gathered, status == "diseased")
dis_list <- split(dis$avg, dis$gene)

##the below is...ugly. But I couldnt' get lapply to work over the two lists. And it was taking too long to figure out how. 
res = NULL
for (i in 1:length(healthy_list)){
  a <- t.test(healthy_list[[i]], dis_list[[i]])$p.value
  res <- c(res, a)
}
##No two genes are sig diff (Welch Two Sample t-test). 
res


#####PLOTS
ggplot(mean_depths_gathered, aes(x = gene, y = avg)) +
  geom_boxplot(aes(fill = status)) +
  theme_classic() +
  scale_fill_grey(start = 0.5) + 
  ylab("Average read depth across all groups within a population") +
  xlab("") +
  theme(axis.text.x = element_text(hjust = 1, angle = 60))
       