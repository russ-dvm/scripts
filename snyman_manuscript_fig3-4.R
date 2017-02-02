##FIGURES: Sex-linked DDX3Y & MBL2 HAPLOTYPES

library(reshape)
library(ggplot2)


#######GENDER PLOT######
gender <- read.table("~/Dropbox/temp/gender.txt", h=T)

##ADD RANKS
gender_sorted <- gender[order(gender["DDX3Y"]),]
gender_sorted$rank <- c(1:nrow(gender_sorted))

##REMOVE OUTLIERS
## PIG 33 == 352; pig 74 == 812
gender_sorted <- subset(gender_sorted, gender_sorted$pig_id != 352)
gender_sorted <- subset(gender_sorted, gender_sorted$pig_id != 812)

##PLOT
ggplot(gender_sorted, aes(x=rank, y=DDX3Y)) + geom_line(colour="light grey") + geom_point(aes(colour = gender), shape = 18, size = 4) + theme_bw() + theme() + xlab("") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line(colour="light grey"), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + scale_y_continuous(breaks=seq(-11,10,1)) + geom_hline(yintercept=0) + ylab(bquote('Gene expression relative to GAPDH ('~log[2]~')')) 

#####MBL2 PLOT######

mbl2_raw <- read.table("~/Dropbox/temp/mbl2_haplotypes.txt", h=T)
mbl2_sorted <- mbl2_raw[order(mbl2_raw["MBL2"]),]
mbl2_sorted$rank <- c(1:nrow(mbl2_sorted))


##REMOVE OUTLIER PIGGIES
mbl2_sorted <- subset(mbl2_sorted, mbl2_sorted$pig_id != 352)
mbl2_sorted <- subset(mbl2_sorted, mbl2_sorted$pig_id != 812)

ggplot(mbl2_sorted, aes(x=rank, y=MBL2)) + geom_line(colour="light grey") + geom_point(aes(colour = haplotype, shape=haplotype), size = 2) + theme_bw() + theme() + xlab("") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line(colour="light grey"), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + scale_y_continuous(breaks=seq(-11,10,1)) + geom_hline(yintercept=0) + ylab(bquote('Gene expression relative to GAPDH ('~log[2]~')')) + scale_shape_manual(values=c(15,16,17,18,19,24))

ggplot(mbl2_sorted, aes(x=rank, y=MBL2)) + 
  geom_line(colour="light grey") + 
  geom_point(aes(shape=haplotype), size = 2) + 
  theme_bw() + 
  xlab("") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line(colour="light grey"), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + scale_y_continuous(breaks=seq(-11,10,1)) + 
  geom_hline(yintercept=0) + 
  ylab(bquote('Gene expression relative to GAPDH ('~log[2]~')')) + 
  scale_shape_manual(values=c(0,1,2,3,4,5))

test <- mbl2_sorted
test$eh <- 1
ggplot(test) + geom_boxplot(aes(x=eh, y=MBL2, fill=haplotype))

