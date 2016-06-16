library(ggplot2)


comparison <- read.table("~/sandbox/genotype_comparison/results_for_R.txt", h=T)

nrow(comparison)
head(comparison)

comp.with.values <- subset(comparison, comparison$maj_af_illumina != "NA")
nrow(comp.with.values)

fully.genotyped <- subset(comp.with.values, comp.with.values$ref_sequenom+comp.with.values$alt_sequenom == 6)

ggplot(fully.genotyped) + geom_jitter(aes(x=maj_af_sequenom, y=maj_af_illumina)) + facet_wrap(~id)
head(fully.genotyped, 100)
##Graphs
ggplot(comp.with.values, aes(x=maj_af_sequenom, y=maj_af_illumina)) + geom_jitter(aes(color=depth_bin)) + geom_smooth(method="lm", se=T) +geom_abline(intercept = 0, slope = 1, color="green") + scale_y_continuous(breaks=seq(0,4,0.25)) + theme_bw() + facet_wrap(~depth_bin)

cor(comp.with.values$maj_af_sequenom, comp.with.values$maj_af_illumina, use="complete")

ggplot(comp.with.values, aes(x=maj_af_sequenom, y=maj_af_illumina)) + geom_jitter(aes(color = depth_bin)) + geom_smooth(method="lm") + facet_wrap(~id) + theme_bw()+geom_abline(intercept = 0, slope = 1, color="green") + scale_y_continuous(breaks=seq(0,4,0.25))


##Diff graphs
genotype.diffs <- comp.with.values

genotype.diffs$diff <- 1 - (abs(test$maj_af_illumina - test$maj_af_sequenom))
genotype.diffs$unique.id <- paste(test$group, test$id, sep = ".")
ggplot(genotype.diffs) + geom_histogram(aes(x=diff, fill=depth_bin), binwidth = (1/6)) + theme_bw()


##without low illumina depth
withmindepth <- subset(comp.with.values, comp.with.values$depth_bin != "low")
ggplot(withmindepth, aes(x=maj_af_sequenom, y=maj_af_illumina)) + geom_jitter(aes(colour=depth_bin)) + geom_smooth(method="lm", se=T) + theme_bw() + geom_abline(intercept = 0, slope = 1, color ="green") + facet_wrap(~id)

cor(withmindepth$maj_af_illumina, withmindepth$maj_af_sequenom)
nrow(withmindepth)
