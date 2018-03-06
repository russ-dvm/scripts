library(tidyverse)
library(gridExtra)
library(viridis)

raw <- read.table("~/porcine/sequenom_results/cleaned_data.txt", h = T, sep = "\t", na.strings = "0", stringsAsFactors = F)
raw$id[1] <- "pig0"

ref <- read.table("~/porcine/sequenom_results/ref_alleles.txt", h = F, sep = "\t")
colnames(ref) <- c("SNP", "ref")


rawTidy <- gather(raw, key = "SNP", value = "geno", colnames(raw)[2:ncol(raw)])

rawTidyRef <- merge(rawTidy, ref)

rawTidyRefSep <- separate(rawTidyRef, geno, c("genoA", "genoB"), sep = 1)
rawTidyRefSep$A <- rawTidyRefSep$genoA != rawTidyRefSep$ref
rawTidyRefSep$B <- rawTidyRefSep$genoB != rawTidyRefSep$ref
rawTidyRefSep$score <- rawTidyRefSep$A + rawTidyRefSep$B
new <- rawTidyRefSep
new$status <- NA
new$status[grep("-", new$id)] <- "diseased"
new$status[is.na(new$status)] <- "normal"


## Frequency of genotypes between diseased and normal pigs - all pigs plotted
ggplot(new, aes(x = score, fill = status, colour = status)) + geom_histogram(aes(y=..density..), position = "dodge", bins = 3, binwidth =1 ,boundary = 0.5) + facet_wrap(~SNP) + xlab("")


## import plink phenotype files to subset data by pathogen
ecoli <- read.table("~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/eschericia_coli_k88.pheno", sep = "\t", h = F, na.strings = "-9")
myco <- read.table("~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/mycoplasma_spp.pheno", sep = "\t", h = F, na.strings = "-9")
flu <- read.table("~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/siv.pheno", sep = "\t", h = F, na.strings = "-9")
prrs <- read.table("~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/prrs.pheno", sep = "\t", h = F, na.strings = "-9")
enteritis <- read.table("~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/enteritis.pheno", sep = "\t", h = F, na.strings = "-9")
sept <- read.table("~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/septicemia.pheno", sep = "\t", h = F, na.strings = "-9")
pneu <- read.table("~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/pneumonia.pheno", sep = "\t", h = F, na.strings = "-9")

## Keep only genotyped animals
ecoli <- ecoli[!is.na(ecoli$V3),]$V1
myco <- myco[!is.na(myco$V3),]$V1
flu <- flu[!is.na(flu$V3),]$V1
prrs <- prrs[!is.na(prrs$V3),]$V1
enteritis <- enteritis[!is.na(enteritis$V3),]$V1
sept <- sept[!is.na(sept$V3),]$V1
pneu <- pneu[!is.na(pneu$V3), ]$V1


## Subset the dataframe that has all the genotypes by animals genotyped for that condition
ecoli_sub <- subset(new, new$id %in% ecoli)
myco_sub <- subset(new, new$id %in% myco)
flu_sub <- subset(new, new$id %in% flu)
prrs_sub <- subset(new, new$id %in% prrs)
enteritis_sub <- subset(new, new$id %in% enteritis)
sept_sub <- subset(new, new$id %in% sept)
pneu_sub <- subset(new, new$id %in% pneu)

## Vectors of SNPs determined to be significant
ecoli_sig <- c("rs330938856","rs322614244","rs327156991","rs330573086","rs319252499","rs335696039","rs343093524","rs327744812","rs343548293","rs346369612","rs339063186","rs336562067")
myco_sig <- c("rs333222079","rs338072079","rs326158227","rs320342887","rs342195622","rs328892780","rs331664783","rs333340587","rs337034769",              "rs330321112","rs322664568","rs335102706","rs81227376")
flu_sig <- c("rs335696039","rs343093524")
prrs_sig <- c("rs339469390","rs333222079","rs338072079")
enteritis_sig <- c("rs335102706","rs330938856","rs322614244")
sept_sig <- c("rs326158227","rs320342887","rs342195622","rs328892780","rs331664783","rs333340587","rs337034769","rs330321112","rs322664568","rs335102706","rs330573086","rs319252499")
pneu_sig <- c("rs326158227","rs320342887","rs342195622","rs331664783","rs337034769","rs330321112","rs322664568","rs335102706","rs330938856","rs322614244","rs330573086","rs319252499","rs343304937","rs320914345")

##Dataframes containing only the a) sig snps and b) animals that were diagnosed with diseases that those snps were significant in
ecoli_sub_sig <- subset(ecoli_sub, ecoli_sub$SNP %in% ecoli_sig)
myco_sub_sig <- subset(myco_sub, myco_sub$SNP %in% myco_sig)
flu_sub_sig <- subset(flu_sub, flu_sub$SNP %in% flu_sig)
prrs_sub_sig <- subset(prrs_sub, prrs_sub$SNP %in% prrs_sig)
enteritis_sub_sig <- subset(enteritis_sub, enteritis_sub$SNP %in% enteritis_sig)
sept_sub_sig <- subset(sept_sub, sept_sub$SNP %in% sept_sig)
pneu_sub_sig <- subset(pneu_sub, pneu_sub$SNP %in% pneu_sig)

##Plots
ecoli.bar <- ggplot(ecoli_sub_sig, aes(x = score, fill = status)) + geom_histogram(aes(y=..density..), position = "dodge", bins = 3, binwidth =1 ,boundary = 0.5) + facet_wrap(~SNP) + theme_bw() + scale_fill_viridis(begin = 0.35, end = 0.8, discrete = T, option = "A") + xlab("") + ylab("Proportion of pigs") + theme(legend.position = "top", legend.title = element_blank())

myco.bar <- ggplot(myco_sub_sig, aes(x = score, fill = status)) + geom_histogram(aes(y=..density..), position = "dodge", bins = 3, binwidth =1 ,boundary = 0.5) + facet_wrap(~SNP) + theme_bw() + scale_fill_viridis(begin = 0.35, end = 0.8, discrete = T, option = "A") + xlab("") + ylab("Proportion of pigs") + theme(legend.position = "top", legend.title = element_blank())

flu.bar <- ggplot(flu_sub_sig, aes(x = score, fill = status)) + geom_histogram(aes(y=..density..), position = "dodge", bins = 3, binwidth =1 ,boundary = 0.5) + facet_wrap(~SNP) + theme_bw() + scale_fill_viridis(begin = 0.35, end = 0.8, discrete = T, option = "A") + xlab("") + ylab("Proportion of pigs") + theme(legend.position = "top", legend.title = element_blank())

prrs.bar <- ggplot(prrs_sub_sig, aes(x = score, fill = status)) + geom_histogram(aes(y=..density..), position = "dodge", bins = 3, binwidth =1 ,boundary = 0.5) + facet_wrap(~SNP) + theme_bw() + scale_fill_viridis(begin = 0.35, end = 0.8, discrete = T, option = "A") + xlab("") + ylab("Proportion of pigs") + theme(legend.position = "top", legend.title = element_blank())

enteritis.bar <- ggplot(enteritis_sub_sig, aes(x = score, fill = status)) + geom_histogram(aes(y=..density..), position = "dodge", bins = 3, binwidth =1 ,boundary = 0.5) + facet_wrap(~SNP) + theme_bw() + scale_fill_viridis(begin = 0.35, end = 0.8, discrete = T, option = "A") + xlab("") + ylab("Proportion of pigs") + theme(legend.position = "top", legend.title = element_blank())

pneu.bar <- ggplot(pneu_sub_sig, aes(x = score, fill = status)) + geom_histogram(aes(y=..density..), position = "dodge", bins = 3, binwidth =1 ,boundary = 0.5) + facet_wrap(~SNP) + theme_bw() + scale_fill_viridis(begin = 0.35, end = 0.8, discrete = T, option = "A") + xlab("") + ylab("Proportion of pigs") + theme(legend.position = "top", legend.title = element_blank())

sept.bar <- ggplot(sept_sub_sig, aes(x = score, fill = status)) + geom_histogram(aes(y=..density..), position = "dodge", bins = 3, binwidth =1 ,boundary = 0.5) + facet_wrap(~SNP) + theme_bw() + scale_fill_viridis(begin = 0.35, end = 0.8, discrete = T, option = "A") + xlab("") + ylab("Proportion of pigs") + theme(legend.position = "top", legend.title = element_blank())


## make eQTL plots to go along with this
#other snp is rs343093524.207.txt
# rs335696039.366.txt
# comparison <- read.table("~/porcine/microarray/eQTL_results_no-outliers/eQTL_results_LR_min3/sequenomed/rs343093524.207.txt", sep = "\t", h = T, nrow = 69)
# gene <- colnames(comparison)
# comparison$snp <- "Rs123"
# gene.df <- as.data.frame(gene)
# gene.name <- gene.df[3,1]
# gene.name <- gsub("expression.", "", gene.name)
# 
# ref <- subset(comparison, genotype == 0, select=genotype_letters)
# ref <- ref[1,]
# hetero <- subset(comparison, genotype == 1, select = genotype_letters)
# hetero <- hetero[1,]
# alt <- subset(comparison, genotype == 2, select = genotype_letters)
# alt <- alt[1,]
# 
# #plot
# plot1 <- ggplot(comparison) + geom_boxplot(aes(x=genotype, y=comparison[,3], group=genotype_letters)) + geom_smooth(aes(x=genotype, y=comparison[,3]), se=F,method="lm") + theme_bw() + geom_jitter(aes(shape=genotype_letters, color=genotype_letters, x=genotype, y=comparison[,3])) + ylab("log2(expression)") + theme(axis.ticks = element_blank()) + theme(legend.position="none") + xlab("") + scale_x_continuous(breaks=seq(0,2,1), labels=c(paste(ref), paste(hetero), paste(alt)))
# plot


## Import all the eQTL data.
file.list <- list.files(path="~/porcine/microarray/eQTL_results_no-outliers/eQTL_results_LR_min3/sequenomed/", pattern=".txt")
path="~/porcine/microarray/eQTL_results_no-outliers/eQTL_results_LR_min3/sequenomed/"

num.files = length(file.list)
mega <- data.frame("animal" = factor(), "genotype" = integer(), "expression" = numeric(), "genotype_letters" = factor())

for (i in c(1:num.files)) {
  
  filename <- paste(path, file.list[i], sep="")
  
    #Import the data into a data frame. Only take the first 72 rows from the file.
  qtl <- read.table(filename, h=T, nrows=69)
  colnames(qtl)[3] <- "expression"
  qtl$snp <- file.list[i]
  mega <- rbind(mega, qtl)

}

## warning will occur because the .number.txt is getting stored.
mega_1 <- separate(mega, snp, c("snp"))
mega_2 <- mega_1[!is.na(mega_1$genotype),]
mega_2$genotype <- as.factor(mega_2$genotype)

box <- ggplot(subset(mega_2, mega_2$snp %in% enteritis_sig), aes(x = genotype, y = expression)) + geom_boxplot(aes(group = genotype_letters)) + geom_jitter(aes(shape=genotype, color=genotype), width = 0.1) + facet_wrap(~snp) + theme_bw() + theme(legend.position = "none") + geom_smooth(method = "lm", se = F) + ylab("log2(expression)")

grid.arrange(enteritis.bar, box)

