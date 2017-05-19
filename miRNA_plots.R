library(ggplot2)
library(viridis)
library(plyr)

##########TARGETSCAN
#miRNA target sites
mirna <- read.table("~/equine/2014_11_24/miRNA/targetscan/equine_miR_analysis_with_coords_and_genes.txt", sep="\t")

mirna$V18 <- factor(mirna$V18, levels = c("2501 - 3000", "2001 - 2500", "1501 - 2000", "1001 - 1500", "501 - 1000", "<= 500"))

#Number of binding sites for each gene
ggplot(mirna, aes(x=V17)) + geom_histogram(stat="count", aes(fill=V18)) + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="RdPu") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of targets") + ggtitle("Number of putative miRNA binding sites") + xlab("")

#Frequency of different miRs
ggplot(mirna, aes(x=V2)) +  geom_histogram(stat="count", aes(fill=V17))+ theme_classic() + scale_fill_brewer(name="Distance from stop codon (bp)", palette="BrBG") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of occurances") + ggtitle("Frequency of different miRs in collectins") + xlab("")

#SNPs in downstream region
all_snps <- read.table("~/equine/2014_11_24/miRNA/targetscan/downstream_3kb_snps_with_genes.txt", sep="\t")
all_snps$V32 <- factor(all_snps$V32, levels = c("2501 - 3000", "2001 - 2500", "1501 - 2000", "1001 - 1500", "501 - 1000", "<= 500"))

ggplot(all_snps, aes(x=V31)) + geom_histogram(aes(fill=V32), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="RdPu") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs in the 3kb downstream from the stop codon") + xlab("")


#SNPs overlapping miRNA targets
snps <- read.table("~/equine/2014_11_24/miRNA/targetscan/snps_in_miR_targets_with_genes.vcf", sep="\t")

snps$V32 <- factor(snps$V32, levels = c("2501 - 3000", "2001 - 2500", "1501 - 2000", "1001 - 1500", "501 - 1000", "<= 500"))

ggplot(snps, aes(x=V31)) + geom_histogram(aes(fill=V32), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="RdPu") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs impacting putative miRNA binding sites (targetscan)") + xlab("")



#####COMBINED#####

default_full <- read.table("~/equine/2014_11_24/miRNA/combined/targetscan_miranda_default_full_with_gene_coords.txt", sep="\t")
default_seed <- read.table("~/equine/2014_11_24/miRNA/combined/targetscan_miranda_default_seed_with_gene_coord.txt", sep="\t")
stringent_full <- read.table("~/equine/2014_11_24/miRNA/combined/targetscan_miranda_stringent_full_with_gene_coord.txt", sep="\t")
stringent_seed <- read.table("~/equine/2014_11_24/miRNA/combined/targetscan_miranda_stringent_seed_with_gene_coord.txt", sep="\t")
compromise_full <- read.table("~/equine/2014_11_24/miRNA/combined/targetscan_miranda_compromise_full_with_gene_coords.txt", sep="\t")
compromise_seed <- read.table("~/equine/2014_11_24/miRNA/combined/targetscan_miranda_compromise_seed_with_gene_coords.txt", sep="\t")

df_list <- list(default_full, default_seed, stringent_full, stringent_seed, compromise_full, compromise_seed)

refactor <- function(x){
  x$V7 <- factor(x$V7, levels = c("2501 - 3000", "2001 - 2500", "1501 - 2000", "1001 - 1500", "501 - 1000", "<= 500"))
  x
}

df_list <- lapply(df_list, refactor)


ggplot(df_list[[1]], aes(x=V6)) + geom_histogram(aes(fill=V7), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="BuGn") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs predicted by both TargetScan and miRanda (default, full miR sequence)") + xlab("")

ggplot(df_list[[2]], aes(x=V6)) + geom_histogram(aes(fill=V7), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="BuGn") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs predicted by both TargetScan and miRanda (default, seed miR sequence)") + xlab("")

ggplot(df_list[[3]], aes(x=V6)) + geom_histogram(aes(fill=V7), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="BuGn") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs predicted by both TargetScan and miRanda (stringent, full miR sequence)") + xlab("")

ggplot(df_list[[4]], aes(x=V6)) + geom_histogram(aes(fill=V7), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="BuGn") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs predicted by both TargetScan and miRanda (stringent, seed miR sequence)") + xlab("")

df_list[[1]]$V8 <- "def_full"
df_list[[2]]$V8 <- "def_seed"
df_list[[3]]$V8 <- "stringent_full"
df_list[[4]]$V8 <- "stringent_seed"
df_list[[5]]$V8 <- "compromise_full"
df_list[[6]]$V8 <- "compromise_seed"

big_frame <- ldply(df_list, rbind)

ggplot(big_frame, aes(x=V6)) + geom_histogram(aes(fill=V8), stat="count", position="dodge") + theme_classic() + theme(legend.position=c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + xlab("") + scale_fill_viridis(discrete=T)

ggplot(big_frame, aes(x=V6)) + geom_histogram(aes(fill=V7), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="BuGn") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs predicted by both TargetScan and miRanda (stringent, seed miR sequence)") + xlab("") + facet_wrap(~ V8)

