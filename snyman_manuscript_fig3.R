#RS Fraser
##FIGURE 3: Sex-linked DDX3Y

library(reshape2)
library(ggplot2)


#######GENDER PLOT######
gender <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/gender.txt", h=T)
gapdh <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/gapdh.txt", h=T)

##MERGE DATAFRAMES AND CLEAN UP -- note, next time use tidyr - easier. 
colnames(gender) <- c("pig_id", "gender", "DDX3Y")
gender$variable <- paste("X", c(1:88), sep="")
gapdh <- melt(gapdh)
gender_merged <- merge(gender, gapdh, by="variable", suffix=c(".GENDER", ".GAPDH"))
gender_merged <- gender_merged[, -5]
gender_merged$normalized <- gender_merged$DDX3Y - gender_merged$value

##ADD RANKS
gender_merged <- gender_merged[order(gender_merged["normalized"]),]
gender_merged$rank <- c(1:nrow(gender_merged))

colnames(gender_merged)[3] <- "Gender"

##REMOVE OUTLIERS
## PIG 33 == 357; pig 74 == 812
gender_merged <- subset(gender_merged, gender_merged$pig_id != 357)
gender_merged <- subset(gender_merged, gender_merged$pig_id != 812)

##PLOT
ggplot(gender_merged, aes(x=rank, y=normalized)) + 
  geom_hline(yintercept=0) + 
  geom_line(colour="light grey") + 
  geom_point(aes(colour = Gender, shape=Gender), size = 2.5) + 
  theme_bw() + 
  xlab("") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line(colour="light grey"), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank(), legend.justification = c(1,0), legend.position = c(0.97,0.65)) + 
  scale_y_continuous(breaks=seq(-11,10,1)) + 
  ylab(expression(paste(italic("DDX3Y "), "expression relative to", italic(" GAPDH "), '('~log[2]~')'))) +
  scale_color_manual(values =c("grey35", "grey71"))
