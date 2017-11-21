library(ggplot2)

snps <- read.table("~/porcine/sequenom_results/genotyping_snp_success.txt", h = T)
samples <- read.table("~/porcine/sequenom_results/genotyping_sample_success.txt", h = T, sep = "\t")

snps$group <- "a"
samples$group <- "a"

min(snps$rate)
min(samples$rate)

summary(snps)
summary(samples)

ggplot(snps, aes(x = group, y = rate)) + geom_boxplot() + theme_bw() + scale_y_continuous(limits = c(0,100))
ggplot(samples, aes(x = group, y = rate)) + geom_boxplot() + theme_bw() + scale_y_continuous(limits = c(0,100)) 
