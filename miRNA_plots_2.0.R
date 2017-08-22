library(ggplot2)

mirna_seed <- read.table("~/Dropbox/temp/miRNA/combined_by_mir_overlap/overlap.snps.named.txt", sep = "\t")
mirna_full <- read.table("~/Dropbox/temp/miRNA/combined_by_mir_overlap/overlap.snps.full.named.txt", sep = "\t")

mirna$V33 <- factor(mirna$V33, levels = c("2501 - 3000", "2001 - 2500", "1501 - 2000", "1001 - 1500", "501 - 1000", "<= 500"))
mirna$V33 <- factor(mirna$V33, levels = c("2501 - 3000", "2001 - 2500", "1501 - 2000", "1001 - 1500", "501 - 1000", "<= 500"))


ggplot(mirna, aes(x=V32)) + 
  geom_histogram(aes(fill = V33), stat = "count") +
  theme_bw() +
  scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="RdPu") + 
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylim(0,10) +
  scale_y_continuous(breaks = seq(0,9,1)) + 
  theme(panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) +
  ylab("Number of variants") +
  ggtitle("Number of variants predicted to fall within the seed region of an miRNA binding element predicted by both TargetScan and miRanda")
       
ggplot(mirna_full, aes(x=V32)) + 
  geom_histogram(aes(fill = V33), stat = "count") +
  theme_bw() +
  scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="RdPu") + 
  xlab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylim(0,10) +
  scale_y_continuous(breaks = seq(0,28,1)) + 
  theme(panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) +
  ylab("Number of variants") +
  ggtitle("Number of variants predicted to fall within an miRNA binding element predicted by both TargetScan and miRanda")
