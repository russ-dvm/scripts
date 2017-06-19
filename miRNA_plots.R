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

default_full <- read.table("~/equine/2014_11_24/miRNA/combined/targetscan_miranda_default_full_final.vcf", sep="\t", h=T)
default_seed <- read.table("~/equine/2014_11_24/miRNA/combined/targetscan_miranda_default_seed_final.vcf", sep="\t", h=T)
stringent_full <- read.table("~/equine/2014_11_24/miRNA/combined/targetscan_miranda_stringent_full_final.vcf", sep="\t", h=T)
stringent_seed <- read.table("~/equine/2014_11_24/miRNA/combined/targetscan_miranda_stringent_seed_final.vcf", sep="\t", h=T)
compromise_full <- read.table("~/equine/2014_11_24/miRNA/combined/targetscan_miranda_compromise_full_final.vcf", sep="\t", h=T)
compromise_seed <- read.table("~/equine/2014_11_24/miRNA/combined/targetscan_miranda_compromise_seed_final.vcf", sep="\t", h=T)

df_list <- list(default_full, default_seed, stringent_full, stringent_seed, compromise_full, compromise_seed)

refactor <- function(x){
  x$DISTANCE <- factor(x$DISTANCE, levels = c("2501 - 3000", "2001 - 2500", "1501 - 2000", "1001 - 1500", "501 - 1000", "<= 500"))
  x
}

df_list <- lapply(df_list, refactor)


ggplot(df_list[[1]], aes(x=GENE)) + geom_histogram(aes(fill=DISTANCE), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="BuGn") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs predicted by both TargetScan and miRanda (default, full miR sequence)") + xlab("")+ scale_y_continuous(breaks=seq(-10,20,1))

ggplot(df_list[[2]], aes(x=GENE)) + geom_histogram(aes(fill=DISTANCE), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="BuGn") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs predicted by both TargetScan and miRanda (default, seed miR sequence)") + xlab("")+ scale_y_continuous(breaks=seq(-10,20,1))

ggplot(df_list[[3]], aes(x=GENE)) + geom_histogram(aes(fill=DISTANCE), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="BuGn") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs predicted by both TargetScan and miRanda (stringent, full miR sequence)") + xlab("")+ scale_y_continuous(breaks=seq(-10,20,1))

ggplot(df_list[[4]], aes(x=GENE)) + geom_histogram(aes(fill=DISTANCE), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="BuGn") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs predicted by both TargetScan and miRanda (stringent, seed miR sequence)") + xlab("")+ scale_y_continuous(breaks=seq(-10,20,1))

ggplot(df_list[[5]], aes(x=GENE)) + geom_histogram(aes(fill=DISTANCE), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="BuGn") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs predicted by both TargetScan and miRanda (compromise, full miR sequence)") + xlab("")+ scale_y_continuous(breaks=seq(-10,20,1))

ggplot(df_list[[6]], aes(x=GENE)) + geom_histogram(aes(fill=DISTANCE), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="BuGn") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs predicted by both TargetScan and miRanda (compromise, seed miR sequence)") + xlab("") + scale_y_continuous(breaks=seq(-10,20,1))


df_list[[1]]$strategy <- "def_full"
df_list[[2]]$strategy <- "def_seed"
df_list[[3]]$strategy <- "stringent_full"
df_list[[4]]$strategy <- "stringent_seed"
df_list[[5]]$strategy <- "compromise_full"
df_list[[6]]$strategy <- "compromise_seed"

big_frame <- ldply(df_list, rbind)

ggplot(big_frame, aes(x=GENE)) + geom_histogram(aes(fill=strategy), stat="count", position="dodge") + theme_classic() + theme(legend.position=c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + xlab("") + scale_fill_viridis(discrete=T) + scale_y_continuous(breaks=seq(-10,20,1))

ggplot(big_frame, aes(x=GENE)) + geom_histogram(aes(fill=DISTANCE), stat="count") + theme_classic() + scale_fill_brewer(breaks=c("<= 500", "501 - 1000", "1001 - 1500", "1501 - 2000", "2001 - 2500", "2501 - 3000"), name="Distance from stop codon (bp)", palette="BuGn") + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("# of SNVs") + ggtitle("Number of SNVs predicted by both TargetScan and miRanda (stringent, seed miR sequence)") + xlab("") + facet_wrap(~ DISTANCE)+ scale_y_continuous(breaks=seq(-10,20,1))

ggplot(big_frame, aes(x=GENE)) + geom_bar(aes(y=miR_count, fill=strategy), stat="identity", position="dodge")



###Tables
aa<-compromise_seed[,c(1:5,31:34)]
write.table(aa, file="~/Dropbox/compromise_seed.txt", sep="\t", row.names = F, quote=F)
            
colnames(compromise_seed)
ncol(compromise_full)

