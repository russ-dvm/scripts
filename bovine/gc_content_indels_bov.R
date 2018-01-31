library(tidyverse)

gc <- read.table("~/bovine/merged_runs/gc_content/gc_content.txt", h = T, stringsAsFactors = F)
gc[nrow(gc)+1,] <- gc[grep("MBL1", gc$X4_usercol),]
gc[12,4] <- "MBL1"
gc[13,4] <- "SFTPA1"

gc <- gc[order(gc$X6_pct_gc, decreasing = T),]


gc$label <- "null"

gc[grep("FCN1", gc$X4_usercol),]$label <- "high"
gc[grep("SFTPA1", gc$X4_usercol),]$label <- "high"
gc[grep("SFTPD", gc$X4_usercol),]$label <- "high"
gc[grep("COLEC12", gc$X4_usercol),]$label <- "low"
gc[grep("MASP2", gc$X4_usercol),]$label <- "low"
gc[grep("COLEC10", gc$X4_usercol),]$label <- "low"
gc[grep("MASP1", gc$X4_usercol),]$label <- "low"

ggplot(subset(gc, gc$label != "null"), aes(x = label, y = X6_pct_gc)) + geom_point() + geom_text(aes(label = X4_usercol), nudge_x = 0.02, hjust = 0)

###**IF rate_of_snps_bov.R** has been run, then you can do the following:

a<-tapply(colecTrimmed$rateKb, colecTrimmed$Gene, summary)
b<-ldply(a)

c <- merge(gc, b, by.x = "X4_usercol", by.y = ".id")

ggplot(c, aes(x = Median, y = X6_pct_gc)) + 
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm", se = F, color = "red") +
  geom_text(aes(label = X4_usercol), hjust = 0, nudge_x = 0.05) +
  xlab("Median variant density (per kb)") +
  ylab("Pct GC Content") +
  annotate("label",label = "paste(italic(R), \" = .36\")", x = 13, y = 0.5, parse = T) +
  coord_cartesian(xlim = c(0,15))

cor.test(c$X6_pct_gc, c$Median)


## indels
vcf <- read.table("~/bovine/merged_runs/6_variants/lectin/INDEL.only.vcf")
vcf <- vcf[,c(1:5)]

vcfGenes <- merge(vcf, gc, by.x = "V1", by.y = "X1_usercol")
vcfGenes <- vcfGenes[,c(1:5,8)]
vcfGenesRates <- merge(vcfGenes, b, by.x = "X4_usercol", by.y = ".id")
vcfGenesRates
d <- as.data.frame(table(vcfGenesRates$X4_usercol))
e <- merge(b, d, by.x = ".id", by.y = "Var1", all = T)
e[is.na(e$Freq),]$Freq <- 0

library(ggrepel)

ggplot(e, aes(x = Median, y = Freq)) + 
  geom_smooth(method = "lm", se = F, color = "red") +
  geom_point() + 
  geom_text_repel(aes(label = .id), hjust = 0, nudge_x = 0.05, check_overlap = T) +
  theme_classic() +
  annotate("label",label = "paste(italic(R), \" = .36\")", x = 13, y = 5, parse = T) +
  xlab("Median variant density (per kb)") +
  ylab("Number of in/dels")
  

cor.test(e$Median, e$Freq)
