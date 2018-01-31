library(gridExtra)
library(tidyverse)

####JUST ONTARGET READS####
##READ IN FILES
group.ontarget <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/2_metrics/on_target/group_on.target.reads.with_lib.txt", h=T)
pig.ontarget <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/2_metrics/on_target/pig_on.target.reads.with_lib.txt", h=T)
group.picard <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/2_metrics/picard/group_summary.txt", h=T)
pig.picard <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/2_metrics/picard/pig_summary.txt", h=T)


##Order ascending
test <- pig.ontarget[order(pig.ontarget$reads),]
test$group <- factor(test$group, levels=test$group)

ordered.group.ontarget <- group.ontarget[order(group.ontarget$reads),]
ordered.group.ontarget$group <- factor(ordered.group.ontarget$group, levels=ordered.group.ontarget$group)
head(ordered.group.ontarget)

##GENERATE GRAPHS
group_bar <- ggplot(ordered.group.ontarget) + geom_bar(aes(x=group, y=reads, fill=lib), stat="identity") + theme_classic() + xlab("") + ggtitle("# of on-target reads, Groups") + theme(axis.text.x=element_text(angle = 90, hjust=1, size=5))
group_bar

wpig_bar <- ggplot(pig.ontarget) + geom_bar(aes(x=group, y=reads, fill=lib), stat="identity") + theme_classic() + ggtitle("# of on-target reads, Pigs")
 ggplot(test) + geom_bar(aes(x=group, y=reads), stat="identity", fill="#164F86") + theme_classic() + ggtitle("# of on-target reads, Pigs") + theme(axis.text.x=element_text(angle = 90, hjust=1, size=5)) + xlab("")
#ggsave( filename="~/Desktop/on_target_reads_ordered_wide.png", w=12.14, h=7.18, units=c("in"))

group_box <- ggplot(group.ontarget) + geom_boxplot(aes(x=lib, y=reads)) + theme_classic() + ggtitle("# of on-target reads, Groups")
pig_box <- ggplot(pig.ontarget) + geom_boxplot(aes(x=lib, y=reads)) + theme_classic() + ggtitle("# of on-target reads, Pigs")

both <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/2_metrics/on_target/both.on-target.txt", h=T)

both_bar <- ggplot(both) + geom_bar(aes(x=group, y=reads, fill=lib), stat="identity") + ggtitle("# of on-target reads, All")
both_box <- ggplot(both) + geom_boxplot(aes(x=lib, y=reads)) + ggtitle("# of on-target reads, All") + theme_classic()

##SAVE GRAPHS
grid.arrange(group_bar, group_box, pig_bar, pig_box, both_bar, both_box, nrow=3, ncol=2)
g <- arrangeGrob(group_bar, group_box, pig_bar, pig_box, both_bar, both_box, nrow=3, ncol=2)
ggsave(file="~/Desktop/on_target_reads.png", g)

####COMPARED TO ALL READS####
group.ontarget$total <- group.picard$TOTAL_READS
pig.ontarget$total <- pig.picard$TOTAL_READS

group.ontarget$percent <- group.ontarget$reads/group.ontarget$total*100
pig.ontarget$percent <- pig.ontarget$reads/pig.ontarget$total*100

group_percent_bar <- ggplot(group.ontarget) + geom_bar(aes(x=group, y=percent, fill=lib), stat="identity", position="dodge") + theme_classic() + ggtitle("Percent reads on-target, groups")
pig_percent_bar <- ggplot(pig.ontarget) + geom_bar(aes(x=group, y=percent, fill=lib), stat="identity") + theme_classic() + ggtitle("Percent reads on-target, pig")

group_percent_boxplot <- ggplot(group.ontarget) + geom_boxplot(aes(x=lib, y=percent)) + theme_classic() + ggtitle("Percent reads on-target, groups")
pig_percent_boxplot <- ggplot(pig.ontarget) + geom_boxplot(aes(x=lib, y=percent)) + theme_classic() + ggtitle("Percent reads on-target, pigs")

grid.arrange(group_percent_bar, group_percent_boxplot, pig_percent_bar, pig_percent_boxplot, nrow=2, ncol=2)
h <- arrangeGrob(group_percent_bar, group_percent_boxplot, pig_percent_bar, pig_percent_boxplot, nrow=2, ncol=2)
ggsave(file="~/Desktop/percent_on_target.png", h)


####Depth of sequencing within targets
depth_raw <- read.table("~/porcine/3_aligned/depth-no-outliers.sample_interval_summary", h=T)
depth <- depth_raw
a <- order(depth$average_coverage)
depth <- depth[a,]
depth$Target <- factor(depth$Target, depth$Target)
depth <- depth[,c(1, grep("mean", colnames(depth)))]
depth_g <- gather(depth, key = "pig", value = "mean", grep("mean", colnames(depth)))

str(depth_g)

ggplot(depth_g, aes(x= Target, y = mean)) + 
  geom_boxplot() + 
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 90, 10)) +
  theme(axis.text.x = element_text(angle = 90)) + 
  geom_hline(yintercept = 10, colour = "red") +
  ylab("Average coverage") +
  xlab("")


depth_anova <- aov(depth_g$mean ~ depth_g$Target)
summary(depth_anova)
TukeyHSD(depth_anova)


#### How does the number of variants called in a region compare to the depth?
library(ggpubr)
library(plyr)
v <- read.table("~/porcine/5_snpeff/variants_with_intervals_appended.txt", h=T)
v_counts <- count(v$interval)
colnames(v_counts) <- c("Target", "freq")
avg_depth <- depth_raw[, c(1,3)]
v_d <- merge(v_counts, avg_depth)
ggplot(v_d, aes(x=average_coverage, y = freq)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson", label.x = 2500, label.y = 2500, hjust = 1, vjust = 0) + 
  ylab("Number of variants") +
  xlab("Average coverage") +
  theme_classic()

## Variants by interval

varInt <- read.table("~/porcine/5_snpeff/variants-by-interval.txt", h = F)
colnames(varInt) <- c("Interval", "Count")
sum(varInt$Count)
varInt <- merge(varInt, avg_depth, by.x = "Interval", by.y = "Target")
varInt <- varInt[order(varInt$average_coverage),]
varInt$Interval <- factor(varInt$Interval, levels = varInt$Interval)
ggplot(varInt, aes(x = Interval, y = Count)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))
