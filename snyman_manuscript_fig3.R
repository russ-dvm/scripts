#RS Fraser
##FIGURE 3: Sex-linked DDX3Y
#####NOTE#####
#The first bit is the original implementation of the figures. Following discussions with BL on 17-07-25, we determined that a) we should be using the same averages as were used in the heatmap figure; and b) the values should NOT be GAPDH normalized. New figure script starts at the end. This is not a huge deal for this figure, since there is only one probe for DDX3Y



library(reshape2)
library(ggplot2)


#######GENDER PLOT######
gender <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/gender.txt", h=T)
# gapdh <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/gapdh.txt", h=T)

##MERGE DATAFRAMES AND CLEAN UP -- note, next time use tidyr - easier. 
colnames(gender) <- c("pig_id", "gender", "DDX3Y")
gender$variable <- paste("X", c(1:88), sep="")
# gapdh <- melt(gapdh)
# gender_merged <- merge(gender, gapdh, by="variable", suffix=c(".GENDER", ".GAPDH"))
# gender_merged <- gender_merged[, -5]
# gender_merged$normalized <- gender_merged$DDX3Y - gender_merged$value

##ADD RANKS
# gender_merged <- gender_merged[order(gender_merged["normalized"]),]
# gender_merged$rank <- c(1:nrow(gender_merged))

gender <- gender[order(gender["DDX3Y"]),]
gender$rank <- c(1:nrow(gender))

# colnames(gender_merged)[3] <- "Gender"
colnames(gender)[2] <- "Gender"
colnames(gender)[3] <- "Value"

##REMOVE OUTLIERS
## PIG 33 == 357; pig 74 == 812
# gender_merged <- subset(gender_merged, gender_merged$pig_id != 357)
# gender_merged <- subset(gender_merged, gender_merged$pig_id != 812)
gender <- subset(gender, pig_id != 357)
gender <- subset(gender, pig_id != 812)

##PLOT
fig3 <- ggplot(gender, aes(x=rank, y=Value)) + 
  geom_hline(yintercept=0) + 
  geom_line() + 
  geom_point(aes(fill = Gender, shape=Gender), size = 2) + 
  scale_shape_manual(values=c(23,21)) +
  theme_bw() + 
  xlab("") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line(colour="light grey"), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank(), legend.justification = c(1,0), legend.position = c(0.97,0.20), legend.text = element_text(size=10), legend.title = element_text(size=10), text = element_text(size=10)) + 
  scale_y_continuous(breaks=seq(-11,10,1)) + 
  ylab(expression(paste(italic("DDX3Y "), "expression by gender ", '('~log[2]~')'))) +
  scale_fill_viridis(discrete = T)
fig3
ggsave(fig3, file="~/Desktop/fig3.eps", width = 84, height = 80, units = "mm")

