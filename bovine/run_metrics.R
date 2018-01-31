library(tidyverse)


run1 <- read.table("~/bovine/2015_12_04/hybrid_selection_metrics.txt", h = T, sep = "\t", na="")
run2 <- read.table("~/bovine/2015_12_13/hybrid_selection_metrics.txt", h = T, sep = "\t", na="")
combined <- read.table("~/bovine/merged_runs/hs_selection/summary.txt", h = T, sep = "\t", na="")

run1_groups <- run1[!is.na(run1$SAMPLE),]
run2_groups <- run2[!is.na(run2$SAMPLE),]
run1_groups$run <- factor(1)
run2_groups$run <- factor(2)

combined$run <- factor("Combined")


both <- rbind(run1_groups, run2_groups)
all <- rbind(both, combined)
bothGathered <- gather(both, key = VARIABLE, value = VALUE, -run, -SAMPLE, -TOTAL_READS)

ggplot(both, aes(x = run, y = TOTAL_READS)) + geom_boxplot()

ggplot(both, aes(x = run, y = PCT_SELECTED_BASES*100)) + 
  geom_boxplot() +
  ylab("Percent reads on target")

ggplot(both, aes(x = run, y = MEAN_TARGET_COVERAGE)) + 
  geom_boxplot() + 
  ylab("Mean target coverage by group")

ggplot(all, aes(x = run, y = MEAN_TARGET_COVERAGE)) + 
  geom_boxplot() +
  ylab("Mean target coverage by group") +
  xlab("Run") +
  theme_bw()

both$SAMPLE <- factor(both$SAMPLE, levels = subset(both, run == "1")[order(subset(both, run == "1")$MEAN_TARGET_COVERAGE, decreasing = T),]$SAMPLE)
ggplot(both, aes(x = SAMPLE, y = MEAN_TARGET_COVERAGE)) + 
  geom_bar(aes(fill = run), stat="identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(size = 0.5, linetype = "solid", colour = "black"))



####
genes <- read.table("~/bovine/merged_runs/depth/refined.intervals.sample_interval_summary", h = T)
## Before dropping overall averages, re-order the factor

genes$gene <- NA
genes[grep("chr11:106773027-106834643", genes$Target),]$gene <- "FCN1"
genes[grep("chr14:47260663-47359061", genes$Target),]$gene <- "COLEC10"
genes[grep("chr16:43449519-43480400", genes$Target),]$gene <- "MASP2"
genes[grep("chr1:80546983-80650824", genes$Target),]$gene <- "MASP1"
genes[grep("chr24:35629347-35866133", genes$Target),]$gene <- "COLEC12"
genes[grep("chr26:6294794-6351278", genes$Target),]$gene <- "MBL2"
genes[grep("chr28:35543221-35605674", genes$Target),]$gene <- "CGN1"
genes[grep("chr28:35619503-35725932", genes$Target),]$gene <- "CL46, CL43"
genes[grep("chr28:35765077-35827384", genes$Target),]$gene <- "SFTPD"
genes[grep("chr28:35837919-35870565", genes$Target),]$gene <- "MBL1, SFTPA1"
genes[grep("chr8:112867045-112896491", genes$Target),]$gene <- "COLEC11"


genes$gene <- factor(genes$gene, levels = genes[order(genes$average_coverage, decreasing = T),]$gene)
genes <- genes[,c(1,grep("mean", colnames(genes)), ncol(genes))]

genesGathered <- gather(genes, key = group,value = mean, -Target, -gene)
genesGathered$group <- gsub("_mean_cvg", "", genesGathered$group)


normal <- paste("group", c(1:8), sep = "")
diseased <- paste("group", c(9:24), sep = "")

## Assign infectious and non-infectious status to the groups
genesGathered$status <- NA

for (i in diseased){
    genesGathered[grep(i, genesGathered$group),]$status <- "Infectious"
    print(i)
}
genesGathered[is.na(genesGathered$status),]$status <- "Non-Infectious"

## split into lists for t testing
healthy <- subset(genesGathered, status == "Non-Infectious")
infec <- subset(genesGathered, status = "Infectious")
healthy_list <- split(healthy$mean, healthy$Target)
dis_list <- split(infec$mean, infec$Target)

## Two sample T test to see whether infectious vs non-infectious coverages are sig different
res <- NULL
for (x in names(healthy_list)) {
  print(x)
  a <- t.test(healthy_list[[x]], dis_list[[x]])$p.value
  res <- c(res, a)
}
res < 0.05

## One way ANOVA or KW to see whether depth of coverage is different between targeted regions
aovResults <- aov(genesGathered$mean ~ genesGathered$gene)
summary(aovResults)
tksd <- TukeyHSD(aovResults)
tksdDF <- as.data.frame(tksd[1])
rownames(tksdDF[tksdDF$genesGathered.gene.p.adj < 0.05,])

kruskal.test(genesGathered$mean ~ genesGathered$Target)
library(dunn.test)
dtRes <- dunn.test(genesGathered$mean, g = genesGathered$gene, list = T, method = "bh")
dtRes$comparisons[dtRes$P.adjusted < 0.05][order(dtRes$comparisons[dtRes$P.adjusted < 0.05])]

###2-way ANOVA testing##
## This is unbalanced (different sample sizes in Infection and Non-infectious) 
library(Rmisc)
ab <- summarySE(genesGathered, measurevar = "mean", groupvars = c("gene", "status"))
ggplot(ab, aes(x = gene, y = mean, colour = status)) + geom_errorbar(aes(ymin=mean-se, ymax = mean+se)) + geom_point(shape = 15, size = 4)
model <- lm(mean ~ status + gene + gene:status, data = genesGathered)
library(car)
Anova(model, type="III", contrasts = c("contr.sum", "contr.poly"))
Anova(model, type="II")

#Check model assumptions
#Residuals shoudl be normal-ish
hist(residuals(model))
plot(fitted(model), residuals(model))

## Print out sig diffs
library(multcompView)
library(lsmeans)
sigLet <- cld(lsmeans(model, "gene", adjust = "tukey"), alpha = 0.05, Letters=letters)
sigLet
sigLet <- as.data.frame(sigLet)
sigLet$.group <- gsub(" ", "", sigLet$.group)
cld(lsmeans(model, "status", adjust = "tukey"), alpha = 0.05, Letters=letters)
sigLet$status <- "Infectious"

#
a<-ggplot(genesGathered, aes(x = gene, y = mean, fill = status)) + 
  geom_boxplot() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(linetype = "solid", colour = "black"), legend.title = element_blank()) +
  scale_fill_manual(values = c("grey78", "grey58")) +
  ylab("Mean depth of coverage") +
  geom_text(data = sigLet, aes(x = gene, y = lsmean*1.3, label = .group)) +
  theme(panel.grid = element_blank()) +
  stat_summary(fun.y = mean, colour = "grey98", geom = "text", label = "....", size = 4, position = position_dodge(width = 0.78)) +
  # stat_summary(fun.y = mean, colour = "black", geom = "point", shape = 18, size = 2.5, position = position_dodge(width = 0.8)) +
  xlab("")
a
ggsave(a, file="~/Dropbox/temp/Fig2.eps", width = 174, height = 90, units = "mm")


#Bar chart is perhaps a better choise as it highlights mean and SE, which is really what the 2-way ANOVA with least-square means post-hoc is testing.
#In order to label the bar chart need to have a status column in the sig let dataframe. arbitrarily chosen between infectious/non-infectious.
# In the end BL, JSL and I agreed to stick with the boxplots, but to add in the mean and be more descriptive in the figure legend

a1 <- ggplot(ab, aes(x = gene, y = mean, fill = status)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=mean-se, ymax = mean+se), position = position_dodge(width = 0.9), width = 0.25) +
  theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(linetype = "solid", colour = "black"), legend.title = element_blank()) +
  scale_fill_manual(values = c("orange", "grey")) +
  ylab("Mean depth of coverage") +
  xlab("") +
  scale_y_continuous(limits = c(0,300), expand = c(0,0)) +
  geom_text(data = sigLet, aes(x = gene, y = lsmean*1.2, label = .group))
a1
ggsave(a1, file="~/Dropbox/temp/fig2.eps", width = 174, height = 90, units = "mm")

b<-ggplot(genesGathered, aes(x = status, y = mean)) + 
  geom_boxplot() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(linetype = "solid", colour = "black"), legend.title = element_blank()) +
  scale_fill_manual(values = c("light grey", "dark grey")) +
  ylab("Mean depth of coverage") + 
  theme(panel.grid = element_blank()) +
  xlab("")
b
library(gridExtra)
grid.arrange(a,b, ncol = 2)             


ggplot(genesGathered, aes(x=mean)) + geom_histogram(bins = 100)
shapiro.test(subset(genesGathered, status == "Non-infectious")$mean)

