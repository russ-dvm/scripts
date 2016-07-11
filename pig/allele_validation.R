library(ggplot2)

alleles <- read.table("~/sandbox/allele_freq_comparison/results_for_R.txt", h=T)

##CALCULATE TOTALS
alleles$total_aff <- alleles$AFF_alt + alleles$AFF_ref
alleles$total_unaff <- alleles$UNAFF_alt + alleles$UNAFF_ref

##CALCULATE ALLELE FREQUENCIES
alleles$ill_MAF_diseased <- alleles$ill_alt_diseased/(alleles$ill_alt_diseased + alleles$ill_ref_diseased)
alleles$ill_MAF_normal <- alleles$ill_alt_normal/(alleles$ill_alt_normal + alleles$ill_ref_normal)

alleles$seq_MAF_diseased <- alleles$seq_alt_diseased/(alleles$seq_alt_diseased + alleles$seq_ref_diseased)
alleles$seq_MAF_normal <- alleles$seq_alt_normal/(alleles$seq_alt_normal + alleles$seq_ref_normal)


###################
##PLOT

ggplot(alleles, aes(x=P_illumina_fisher, y=P_seq_fisher)) + geom_point() + geom_smooth() + theme_bw()


ggplot(alleles,aes(x=seq_MAF_diseased, y=ill_MAF_diseased)) + geom_point() + theme_bw() + geom_smooth(method="lm", se=T) +geom_abline(intercept = 0, slope = 1, color="green")


cor(alleles$seq_MAF_diseased, alleles$ill_MAF_diseased, use="complete")

alleles$total_ill_normal_alleles <- alleles$ill_alt_normal + alleles$ill_ref_normal 

deeply.sequenced <- subset(alleles, alleles$total_ill_normal_alleles > 240)

cor(deeply.sequenced$ill_MAF_normal, deeply.sequenced$seq_MAF_normal)
cor(deeply.sequenced$ill_MAF_diseased, deeply.sequenced$seq_MAF_diseased)
cor(deeply.sequenced$P_seq_chisq, deeply.sequenced$P_illumina_chisq)

ggplot(deeply.sequenced, aes(x=seq_MAF_normal, y=ill_MAF_normal)) + geom_point() + theme_bw() + geom_smooth(method="lm", se=T) +geom_abline(intercept = 0, slope = 1, color="green") + geom_point(aes(x=P_illumina_fisher, y=P_seq_fisher), color = "red")

ggplot(deeply.sequenced, aes(x=P_illumina_fisher, y=P_seq_fisher)) + geom_point() + theme_bw() + geom_smooth(method="lm", se=T) + geom_abline(intercept = 0, slope =1, color = "green") + geom_text(aes(label=SNP))


ggplot(alleles, aes(x=seq_MAF_diseased, y=ill_MAF_diseased)) + geom_point()



#####POPOOLATION

ggplot(deeply.sequenced, aes(x=P_seq_chisq, y=P_popoolation)) + geom_point()
head(alleles)
