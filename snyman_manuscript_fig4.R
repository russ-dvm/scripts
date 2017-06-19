#RS Fraser

library(ggplot2)
library(reshape2)
library(viridis)

#####MBL2 PLOT######
gapdh <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/gapdh.txt", h=T)
gapdh <- melt(gapdh)
mbl2_raw <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/mbl2_haplotypes.txt", h=T)
mbl2_raw$variable <- paste("X", c(1:88), sep="")
mbl2_merged <- merge(mbl2_raw, gapdh, by = "variable", suffix = c(".MBL2", ".GAPDH"))
mbl2_merged <- mbl2_merged[, -7]
mbl2_merged$normalized <- mbl2_merged$MBL2 - mbl2_merged$value
mbl2_merged <- mbl2_merged[order(mbl2_merged["normalized"]),]
mbl2_merged$rank <- c(1:nrow(mbl2_merged))
#colnames(mbl2_merged) <- c("Pig", "Gene", "orig_value", "GAPDH", "normalized", "rank")


##REMOVE OUTLIER PIGGIES
mbl2_merged <- subset(mbl2_merged, mbl2_merged$pig_id != 357)
mbl2_merged <- subset(mbl2_merged, mbl2_merged$pig_id != 812)

colnames(mbl2_merged)[5] <- "Haplotype"


###WITH COLOUR
ggplot(mbl2_merged, aes(x=rank, y=normalized)) + 
  geom_line(colour="light grey") + 
  geom_point(aes(colour = Haplotype), size = 2.5) + 
  theme_bw() + theme() +
  xlab("") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line(colour="light grey"), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + scale_y_continuous(breaks=seq(-11,10,1)) + 
  geom_hline(yintercept=0) + 
  ylab(expression(paste(italic("MBL2 "), "expression relative to", italic(" GAPDH "), '('~log[2]~')'))) + 
  scale_shape_manual(values=c(15,16,17,18,19,24)) + 
  scale_color_viridis(discrete = T) 
  ggtitle("MBL2") + 
  theme(plot.title = element_text(hjust = 0.5))

####BY SHAPE
ggplot(mbl2_merged, aes(x=rank, y=normalized)) + 
  geom_line(colour="light grey") + 
  geom_point(aes(shape=Haplotype), size = 2) + 
  theme_bw() + 
  xlab("") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line(colour="light grey"), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + scale_y_continuous(breaks=seq(-11,10,1)) + 
  geom_hline(yintercept=0) + 
  ylab(expression(paste(italic("MBL2 "), "expression relative to", italic(" GAPDH "), '('~log[2]~')'))) + 
  scale_shape_manual(values=c(0,1,2,3,4,5))


