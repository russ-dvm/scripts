library(tidyverse)

raw <- read.table("~/porcine/sequenom_results/cleaned_data.txt", h = T, sep = "\t", na.strings = "0", stringsAsFactors = F)
raw$id[1] <- "pig0"

ref <- read.table("~/porcine/sequenom_results/ref_alleles.txt", h = F, sep = "\t")
colnames(ref) <- c("SNP", "ref")


rawTidy <- gather(raw, key = "SNP", value = "geno", colnames(raw)[2:ncol(raw)])

rawTidyRef <- merge(rawTidy, ref)

rawTidyRefSep <- separate(rawTidyRef, geno, c("genoA", "genoB"), sep = 1)
rawTidyRefSep$A <- rawTidyRefSep$genoA == rawTidyRefSep$ref
rawTidyRefSep$B <- rawTidyRefSep$genoB == rawTidyRefSep$ref
rawTidyRefSep$score <- rawTidyRefSep$A + rawTidyRefSep$B
head(rawTidyRefSep)
new <- rawTidyRefSep
new$status <- NA
new$status[grep("-", new$id)] <- "diseased"
new$status[is.na(new$status)] <- "normal"


## Frequency of genotypes between diseased and normal pigs
ggplot(new, aes(x = score, fill = status, colour = status)) + geom_histogram(aes(y=..density..), position = "dodge", bins = 3, binwidth =1 ,boundary = 0.5) + facet_wrap(~SNP) + xlab("")


## import plink phenotype files to subset data by pathogen
ecoli <- read.table("~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/eschericia_coli_k88.pheno", sep = "\t", h = F, na.strings = "-9")
myco <- read.table("~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/mycoplasma_spp.pheno", sep = "\t", h = F, na.strings = "-9")
flu <- read.table("~/porcine/sequenom_results/plink_analysis/by_condition_or_pathogen/siv.pheno", sep = "\t", h = F, na.strings = "-9")

ecoli <- ecoli[is.na(ecoli$V3) == F,]$V1
myco <- myco[is.na(myco$V3) == F,]$V1
flu <- flu[is.na(flu$V3) == F,]$V1

ecoli_sub <- subset(new, new$id %in% ecoli)
myco_sub <- subset(new, new$id %in% myco)
flu_sub <- subset(new, new$id %in% flu)

ecoli_sig <- c("rs330938856","rs322614244","rs327156991","rs330573086","rs319252499","rs335696039","rs343093524","rs327744812","rs343548293","rs346369612","rs339063186","rs336562067")
myco_sig <- c("rs333222079","rs338072079","rs326158227","rs320342887","rs342195622","rs328892780","rs331664783","rs333340587","rs337034769",              "rs330321112","rs322664568","rs335102706","rs81227376")
flu_sig <- c("rs335696039","rs343093524")

ecoli_sub_sig <- subset(ecoli_sub, ecoli_sub$SNP %in% ecoli_sig)
myco_sub_sig <- subset(myco_sub, myco_sub$SNP %in% myco_sig)
flu_sub_sig <- subset(flu_sub, flu_sub$SNP %in% flu_sig)



ggplot(ecoli_sub_sig, aes(x = score, fill = status, colour = status)) + geom_histogram(aes(y=..density..), position = "dodge", bins = 3, binwidth =1 ,boundary = 0.5) + facet_wrap(~SNP)

ggplot(myco_sub_sig, aes(x = score, fill = status, colour = status)) + geom_histogram(aes(y=..density..), position = "dodge", bins = 3, binwidth =1 ,boundary = 0.5) + facet_wrap(~SNP)

ggplot(flu_sub_sig, aes(x = score, fill = status, colour = status)) + geom_histogram(aes(y=..density..), position = "dodge", bins = 3, binwidth =1 ,boundary = 0.5) + facet_wrap(~SNP)

## make eQTL plots to go along with this
#other snp is rs343093524.207.txt
# rs335696039.366.txt
comparison <- read.table("~/porcine/microarray/eQTL_results_no-outliers/eQTL_results_LR_min3/rs335696039.366.txt", sep = "\t", h = T, nrow = 69)
gene <- colnames(comparison)
gene.df <- as.data.frame(gene)
gene.name <- gene.df[3,1]
gene.name <- gsub("expression.", "", gene.name)

ref <- subset(comparison, genotype == 0, select=genotype_letters)
ref <- ref[1,]
hetero <- subset(comparison, genotype == 1, select = genotype_letters)
hetero <- hetero[1,]
alt <- subset(comparison, genotype == 2, select = genotype_letters)
alt <- alt[1,]

#number of animals genotyped
num.genotyped <- nrow(na.omit(comparison))
annotation <- grobTree(textGrob(paste("# genotyped:", num.genotyped, sep=" "), x=0.9, y=0.98))


#plot
plot <- ggplot(comparison) + geom_boxplot(aes(x=genotype, y=comparison[,3], group=genotype_letters)) + geom_smooth(aes(x=genotype, y=comparison[,3]), se=F,method="lm") + theme_bw() + geom_jitter(aes(shape=genotype_letters, color=genotype_letters, x=genotype, y=comparison[,3])) + ylab("log2(expression)") + ggtitle(paste("Genotype vs Expression for", gene.name)) + theme(axis.ticks = element_blank()) + theme(legend.position="none") + xlab("") + scale_x_continuous(breaks=seq(0,2,1), labels=c(paste(ref), paste(hetero), paste(alt))) + annotation_custom(annotation)
plot
