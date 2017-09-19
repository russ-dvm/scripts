library(ggplot2)
library(ggpubr)

#### IMPORT DATA ####
mirna_seed <- read.table("~/equine/2014_11_24/miRNA/combined_by_mir_overlap/overlap.snps.seed.named.txt", sep = "\t")
mirna_full <- read.table("~/equine/2014_11_24/miRNA/combined_by_mir_overlap/overlap.snps.full.named.txt", sep = "\t")
mirna_ann <- read.table("~/equine/2014_11_24/miRNA/combined_by_mir_overlap/overlap.seed_or_full.txt", sep = "\t")
mirna_sites <- read.table("~/equine/2014_11_24/miRNA/combined_by_mir_overlap/overlap.full.bed", sep = "\t")

mirna_seed$V33 <- factor(mirna_seed$V33, levels = c("2501 - 3000", "2001 - 2500", "1501 - 2000", "1001 - 1500", "501 - 1000", "<= 500"))
mirna_full$V33 <- factor(mirna_full$V33, levels = c("2501 - 3000", "2001 - 2500", "1501 - 2000", "1001 - 1500", "501 - 1000", "<= 500"))
mirna_ann$V33 <- factor(mirna_ann$V33, levels = c("2501 - 3000", "2001 - 2500", "1501 - 2000", "1001 - 1500", "501 - 1000", "<= 500"))


#### PLOTS ####
## SNPS THAT FALL INTO ONLY THE SEED AREA OF THE MRE
ggplot(mirna_seed, aes(x=V32)) + 
  geom_histogram(aes(fill = V33), stat = "count") +
  theme_bw() +
  # scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="RdPu") + 
  scale_fill_viridis(discrete = T, breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", direction = -1) +
  
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylim(0,10) +
  scale_y_continuous(breaks = seq(0,9,1)) + 
  theme(panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) +
  ylab("Number of variants") +
  ggtitle("Number of variants predicted to fall within the seed region of an miRNA binding element predicted by both TargetScan and miRanda")
      
## SNPS THAT FALL INTO THE FULL RANGE OF THE MRE
ggplot(mirna_full, aes(x=V32)) + 
  geom_histogram(aes(fill = V33), stat = "count") +
  theme_bw() +
  # scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="RdPu") + 
  scale_fill_viridis(discrete = T, breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", direction = -1) +
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_y_continuous(breaks = seq(0,28,1)) + 
  theme(panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) +
  ylab("Number of variants") +
  ggtitle("Number of variants predicted to fall within an miRNA binding element predicted by both TargetScan and miRanda") + 
  theme(legend.background = element_rect(fill = "white", colour = "black", size = 0.3, linetype = "solid")) +
  theme(legend.position = c(1,1), legend.justification= c(1,1))


###SEED OR FULL
simple <- data.frame("type" = factor(mirna_ann$V34), "Gene" = factor(mirna_ann$V32), "mirna" = factor(mirna_ann$V31))
missing <- data.frame("type" = NA, "Gene" = c("COLEC12", "MASP1", "MASP2"), "mirna" = NA)
simple <- rbind(simple, missing)

simple$Gene <- factor(simple$Gene, levels = c("COLEC10", "COLEC11", "COLEC12", "FCN1", "FCN1-like", "FCN3", "MASP1", "MASP2", "MBL1", "MBL2", "SFTPA1", "SFTPD"))

ggplot(simple, aes(x = Gene)) + 
  geom_histogram(aes(fill = type), stat = "count", position = "dodge") +
  theme_bw() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylim(0,10) +
  scale_y_continuous(breaks = seq(0,28,1)) + 
  theme(panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) +
  scale_fill_grey(start = 0.5, na.value = "white",
                     name = "Location of variant\nrelative to MRE", 
                     labels = c("Full", "Seed"), 
                     breaks = c("full", "seed")) +
  ylab("Number of variants") +
  ggtitle("") +
  theme(legend.position = c(1,1), legend.justification= c(1,1)) +
  theme(legend.background = element_rect(fill = "white", linetype = "solid", colour = "black", size = 0.25))
    

#### COMP TO TOTAL NUM OF SNPS
## Run "rate_of_snps.R" to get the colec table
dc <- colec[Region == "Downstream 3 kb"]

## Obtain number of SNPs hitting a miRNA seed site for each gene
sm <- count(mirna_seed$V32)
colnames(sm) <- c("Gene", "total")

## Merge the DFs on Gene
q <- as.data.table(merge(dc, sm, all = T))
q[is.na(total)]$total <- 0

rval1 <- cor.test(q$total,q$No.variants)

ggplot(q, aes(x= total, y = No.variants)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  stat_cor(label.x = 2.5)
  ggtitle("SNPs hitting a seed region vs total SNPs")


## Use the annotation_info df from 'rate_of_snps.R' to annotate the mirna_sites file
annotation_info <- as.data.table(annotation_info)
ai <- annotation_info[feature == "downstream_3"]

mirna_sites$avg <- (mirna_sites$V3 + mirna_sites$V2)/2

for (i in 1:nrow(mirna_sites)) {
  ab <- mirna_sites[i,6] >= ai$feature_start
  ba <- mirna_sites[i,6] <= ai$feature_end
  mirna_sites$gene[i] <- ai[ab == ba, ]$gene_name
}

mis <- count(mirna_sites$gene)
colnames(mis) <- c("Gene", "total.mre")
tu <- as.data.table(merge(mis, sm, all = T))
tu[is.na(total)]$total <- 0


ggplot(tu, aes(x=total.mre, y = total)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(label.x = 30) + 
  ggtitle("Total MREs vs SNPs in those MREs")
