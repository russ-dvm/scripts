library(tidyverse)

tf <- read.table("~/bovine/merged_runs/tfbs/intersection_meme-cis/all.summarized", h = F, sep = "\t")

## Do the number of variants in conserved motifs correlate to the total number of variants found in 5kb of that gene?
## Make sure "rate_of_snps_bov.R" has been run first
tf_count <- as.data.frame(table(tf$V24))

combined <- merge(tf_count, subset(colecTrimmedByGene, Region == "5kb_from_start"))
ggplot(combined, aes(x = FReq, y = No.variants)) + geom_point() + geom_smooth(method = "lm")
cor(combined$Freq, combined$No.variants)

tf_levels <- as.data.frame(table(tf$V12))
tf$V12 <- factor(tf$V12, levels = tf_levels[order(tf_levels$Freq, decreasing = T),]$Var1)

tf$V24 <- factor(tf$V24, levels = tf_count[order(tf_count$Freq, decreasing = T), ]$Var1)

## Check the nebulously named "Homeodomain" family
homeo <- subset(tf, tf$V12 == "Homeodomain")


a <- ggplot(tf, aes(x = V24)) + 
  geom_histogram(stat="count") + 
  theme_classic() + 
  theme(text = element_text(size = 11)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab("") +
  ylab("Number of predicted TFBS") +
  scale_y_continuous(expand = c(0,0), limits = c(0,30))
a

b <- ggplot(tf, aes(x = V12)) + 
  geom_histogram(stat="count") + 
  theme_classic() + 
  theme(text = element_text(size = 11)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab("") +
  ylab("Number of occurances") +
  scale_y_continuous(expand = c(0,0), limits = c(0,40))
b

library(gridExtra)
grid.arrange(a, b, nrow = 2)

all <- read.table("~/bovine/merged_runs/tfbs/cis-bp/all.condensed", h = F, sep = "\t")
ggplot(all, aes(x = V19)) + geom_histogram(stat = "count")


## Looking at individual transcription factors
library(data.table)
tfId <- fread("~/bovine/merged_runs/tfbs/cis-bp/all.forward-only.csv", quote = "\"")
test2[order(test2$Freq, decreasing = T), ]