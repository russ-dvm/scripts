## NEED TO RUN rate_of_snps_bov.R and rate_of_snps.R FIRST!!!

cow <- colecTrimmed
horse <- colec
horse$rateKb <- horse$Rate*1000

horse$Species <- "horse"
cow$Species <- "cattle"

combo <- rbind(horse, cow)
combo <- combo[order(combo$rateKb, decreasing = T),]
combo$Gene <- factor(combo$Gene, levels = c("FCN1-like", "FCN1", "SFTPD", "SFTPA1", "CL46", "MBL1", "MBL2", "COLEC11", "CL43", "CGN1", "MASP2","COLEC10",  "FCN3", "COLEC12", "MASP1"))

ggplot(combo, aes(x = Gene, y = rateKb)) + geom_boxplot(aes(fill = Species)) + theme_classic() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + scale_fill_manual(values = c("palegreen3", "palegoldenrod")) + theme(legend.position = c(1,1), legend.justification = c(1,1)) + ylab("Variant density (per kb)") + xlab("")

