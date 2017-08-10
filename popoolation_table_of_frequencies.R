library(tidyverse)

freq_table <- read.table("~/equine/2014_11_24/popoolation/subsampled340_rc", h=T, stringsAsFactors = F)
freq_table_2 <- read.table("~/equine/2014_11_24/popoolation/subsampled340-min20-noindels-for-rc_rc", h=T, stringsAsFactors = F)
fisher <- read.table("~/equine/2014_11_24/popoolation/subsampled340-min20.sig-only.fet", col.names=c("Chrom", "pos", "A", "B", "C", "Pop", "invlogp"))

#Check the the major and minor alleles for each population are the same
freq_table <- separate(freq_table, major_alleles.maa., c("maa_n", "maa_d"), sep = 1)
freq_table <- separate(freq_table, minor_alleles.mia., c("mia_n", "mia_d"), sep = 1)
freq_table$maa_same <- freq_table$maa_n == freq_table$maa_d
freq_table$mia_same <- freq_table$mia_n == freq_table$mia_d
temp <- freq_table

##Fix table where major/minor alleles are different for the two pops
for (i in 1:nrow(freq_table)){
  if (freq_table[i,]$maa_same == F){
    storage_maa_d <- freq_table[i,]$maa_d
    storage_maa_2 <- freq_table[i,]$maa_2
    freq_table[i,]$maa_d <- freq_table[i,]$mia_d
    freq_table[i,]$maa_2 <- freq_table[i,]$mia_2
    freq_table[i,]$mia_d <- storage_maa_d
    freq_table[i,]$mia_2 <- storage_maa_2
  }
  
}

#Re-check that the values are the same
freq_table$maa_same <- freq_table$maa_n == freq_table$maa_d
freq_table$mia_same <- freq_table$mia_n == freq_table$mia_d

##Good idea to inspect the table manually -- for the equine data there are a few SNPs that should be excluded and a couple that need to be manually fixed.
#Entry 542 - normal population was entirely reference, thus the minor allele was not populated (weird)
freq_table[542,]$mia_n <- "A"
freq_table[542,]$mia_d <- "A"
freq_table[2644,]$mia_d <- "A"
freq_table$maa_same <- freq_table$maa_n == freq_table$maa_d
freq_table$mia_same <- freq_table$mia_n == freq_table$mia_d
freq_table <- freq_table[freq_table$mia_same == T,]

##The above should remove these three
# freq_table[-2557,]
# freq_table[-3937,]
# freq_table[-5556,]

##Merge the results of fisher's exact test with the actual frequencies
mega_table <- merge(freq_table, fisher, by.x = "pos")
mega_table_curated <- mega_table[c(-3, -4, -6, -7, -16:-22)]


##Play with frequency diffs
maybe <- data.frame("mia_n" = as.character(mega_table$mia_1), "mia_d" = as.character(mega_table$mia_2))
maybe <- separate(maybe, mia_n, c("alt_n", "total_n"), sep = "/")
maybe <- separate(maybe, mia_d, c("alt_d", "total_d"), sep = "/")
maybe$alt_n <- as.numeric(maybe$alt_n)
maybe$alt_d <- as.numeric(maybe$alt_d)
maybe$total_n <- as.numeric(maybe$total_n)
maybe$total_d <- as.numeric(maybe$total_d)

maybe$n_freq <- maybe$alt_n/maybe$total_n
maybe$d_freq <- maybe$alt_d/maybe$total_d

maybe$diff <- maybe$n_freq - maybe$d_freq
library(MASS)
ggplot(maybe, aes(x = diff)) + geom_histogram(colour = "black", fill = "white") + theme_bw()
qqnorm(maybe$diff) + qqline(maybe$diff)
?qqnorm
