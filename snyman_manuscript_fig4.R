#RS Fraser
#####NOTE#####
#The first bit is the original implementation of the figures. Following discussions with BL on 17-07-25, we determined that a) we should be using the same averages as were used in the heatmap figure; and b) the values should NOT be GAPDH normalized. New figure script starts at the end. This is not a huge deal for this figure, since there is only one probe for MBL2


library(ggplot2)
library(reshape2)
library(viridis)

#####MBL2 PLOT######
# gapdh <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/gapdh.txt", h=T)
# gapdh <- melt(gapdh)
mbl2_raw <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/mbl2_haplotypes.txt", h=T)
# mbl2_raw$variable <- paste("X", c(1:88), sep="")
# mbl2_merged <- merge(mbl2_raw, gapdh, by = "variable", suffix = c(".MBL2", ".GAPDH"))
# mbl2_merged <- mbl2_merged[, -7]
# mbl2_merged$normalized <- mbl2_merged$MBL2 - mbl2_merged$value
# mbl2_merged <- mbl2_merged[order(mbl2_merged["normalized"]),]
mbl2_raw <- mbl2_raw[order(mbl2_raw["MBL2"]),]
mbl2_raw$rank <- c(1:nrow(mbl2_raw))
# mbl2_merged$rank <- c(1:nrow(mbl2_merged))
## colnames(mbl2_merged) <- c("Pig", "Gene", "orig_value", "GAPDH", "normalized", "rank")


##REMOVE OUTLIER PIGGIES
# mbl2_merged <- subset(mbl2_merged, mbl2_merged$pig_id != 357)
# mbl2_merged <- subset(mbl2_merged, mbl2_merged$pig_id != 812)
mbl2_raw <- subset(mbl2_raw, pig_id != 357)
mbl2_raw <- subset(mbl2_raw, pig_id != 812)

# colnames(mbl2_merged)[5] <- "Haplotype"
colnames(mbl2_raw)[4] <- "Haplotype"


###WITH COLOUR
fig4 <- ggplot(mbl2_raw, aes(x=rank, y=MBL2)) + 
  geom_hline(yintercept=0) + 
  geom_line() + 
  geom_point(aes(fill = Haplotype, shape = Haplotype), size = 2) +
  scale_shape_manual(values=c(22,24,21,23,25,22)) +
  theme_bw() +
  xlab("") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line(colour="light grey"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank(), text = element_text(size=10)) + 
  scale_y_continuous(breaks=seq(-11,10,1)) + 
  ylab(expression(paste(italic("MBL2 "), "expression by haplotype ", '('~log[2]~')'))) + 
  scale_fill_viridis(discrete = T) +
  theme(legend.justification = c(1,0), legend.position = c(0.97,0.01), legend.text=element_text(size=10))
fig4
ggsave(fig4, file="~/Desktop/fig4.eps", width = 174, height = 80, units = "mm")

####BY SHAPE
# ggplot(mbl2_merged, aes(x=rank, y=normalized)) + 
#   geom_line(colour="light grey") + 
#   geom_point(aes(shape=Haplotype), size = 2) + 
#   theme_bw() + 
#   xlab("") + 
#   theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line(colour="light grey"), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + scale_y_continuous(breaks=seq(-11,10,1)) + 
#   geom_hline(yintercept=0) + 
#   ylab(expression(paste(italic("MBL2 "), "expression relative to", italic(" GAPDH "), '('~log[2]~')'))) + 
#   scale_shape_manual(values=c(0,1,2,3,4,5))


