library(ggplot2)

predicts <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/3_variants/coding_snp_predictions/polyphen_predictions_individual_pigs.txt", h=T)

predicts <- subset(predicts, predicts$pph2_prob != "?")
predicts <- subset(predicts, predicts$oneminussift != "N/A")

ggplot(predicts, aes(x=pph2_prob, y=oneminussift)) + geom_point() + theme_bw() + theme(axis.text.x = element_text(angle=90, hjust=1, size=5), axis.text.y = element_text(size=6))


str(predicts)
predicts$sift <- as.integer(predicts$SIFT_score)
head(predicts)
