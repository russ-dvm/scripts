library(reshape)
library(ggplot2)

###IMPORT REF GENE (GAPDH)
gapdh <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/gapdh.txt", h=T)
gapdh <- melt(gapdh)


###IMPORT HAPTOGLOBIN -- multiple probes. Disentangle.
x = read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/haptoglobin.txt", h = T)

y <- melt(x)

y <- y[order(y["value"], decreasing = T), ]
z <- y
z$rank <- c(1:nrow(y))
#ggplot(z, aes(x = variable, y = value)) + geom_boxplot() + theme_classic() + xlab("") 
z
#subset based on probe
p1 <- subset(y, y$probe == "A_72_P178006")
p2 <- subset(y, y$probe == "A_72_P568194")
p3 <- subset(y, y$probe == "A_72_P581127")
p4 <- subset(y, y$probe == "A_72_P657000")
p5 <- subset(y, y$probe == "A_72_P667554")
  
#define "order" based on the rank of P1.
p1_sorted <- p1[order(p1["value"]),]
p1_sorted$rank <- c(1:nrow(p1_sorted))
p1_sorted
p1 <- p1_sorted[order(p1_sorted["variable"]),]
default_rank <- p1$rank

#sort by pig # and then give rank
p2 <- p2[order(p2["variable"]),]
p2$rank <- default_rank

p3 <- p3[order(p3["variable"]),]
p3$rank <- default_rank

p4 <- p4[order(p4["variable"]),]
p4$rank <- default_rank

p5 <- p5[order(p5["variable"]),]
p5$rank <- default_rank


complete <- rbind(p1, p2, p3, p4, p5)

#all data for haptoglobin
#ggplot(complete) + geom_point(aes(x=rank, y=value)) + facet_grid(.~probe) + ggtitle("all probes, all pigs")

#drop the lowest two
#ggplot(subset(complete, rank != 1 & rank != 2)) + geom_point(aes(x=rank, y=value)) + facet_grid(.~probe) + ggtitle("all probes, bottom 2 pigs removed")

#average the five probes
avg <- data.frame(Pig = as.character(), value=as.integer())
total = 0
for (i in 1:nrow(p1)){
  rank <- subset(complete, rank==i)
  avg <- rbind(avg, c(rank$variable[1], as.numeric(mean(rank$value))))
  }

#avg <- avg[order(avg[,2]),]
#avg$rank <- c(1:nrow(avg))

#ggplot(avg, aes(x=rank, y=avg[,2])) + geom_point() + ggtitle("average of all 5 probes, all pigs")

#avg_no_out <- subset(avg, rank != 1 & rank != 2)
#ggplot(avg_no_out) + geom_point(aes(x=rank, y=avg_no_out[,2]), shape=5)+ ggtitle("average of all 5 probes, bottom 2 pigs removed")

###2x2
#pub <- rbind(p1, p2, p3, p4)
#ggplot(pub, aes(x=rank, y=value)) + geom_point() + facet_wrap(~probe, scales="free") + theme_bw() + ylab("GAPDH") + xlab("Pig ID")

##
avg$Gene <- "HP"
colnames(avg) <- c("variable", "value", "Gene")
hp_sorted <- data.frame(Gene=avg$Gene, variable = avg$variable, value = avg$value)
hp_sorted$variable = paste("X", hp_sorted$variable, sep = "")
hp_merged <- merge(hp_sorted, gapdh, by = "variable", suffixes = c(".HP", ".GAPDH"))
hp_merged <- hp_merged[, -4]
hp_merged$normalized <- hp_merged$value.HP - hp_merged$value.GAPDH

hp_merged <- hp_merged[order(hp_merged$normalized),]
hp_merged$rank <- c(1:nrow(hp_merged))
colnames(hp_merged) <- c("Pig", "Gene", "orig_value", "GAPDH", "normalized", "rank")


####IMPORT MBL2 & NLRP5

mbl2_raw <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/mbl2.txt", h=T)
mbl2_melted <- melt(mbl2_raw)
mbl2_merged <- merge(mbl2_melted, gapdh, by = "variable", suffix = c(".MBL2", ".GAPDH"))
mbl2_merged <- mbl2_merged[, -4]
mbl2_merged$normalized <- mbl2_merged$value.MBL2 - mbl2_merged$value.GAPDH
mbl2_merged <- mbl2_merged[order(mbl2_merged["normalized"]),]
mbl2_merged$rank <- c(1:nrow(mbl2_merged))
colnames(mbl2_merged) <- c("Pig", "Gene", "orig_value", "GAPDH", "normalized", "rank")


nlrp5 <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/nlrp5.txt", h=T)
nlrp5_melted <- melt(nlrp5)
nlrp5_merged <- merge(nlrp5_melted, gapdh, by="variable", suffix=c(".NLRP5", ".GAPDH"))
nlrp5_merged <- nlrp5_merged[,-4]
nlrp5_merged$normalized <- nlrp5_merged$value.NLRP5 - nlrp5_merged$value.GAPDH
nlrp5_merged <- nlrp5_merged[order(nlrp5_merged["normalized"]),]
nlrp5_merged$rank <- c(1:nrow(nlrp5_merged))
colnames(nlrp5_merged) <- c("Pig", "Gene", "orig_value", "GAPDH", "normalized", "rank")

##HAMP

hamp_raw <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/hamp.txt", h=T)
hamp_melted <- melt(hamp_raw)
hamp_avg <- data.frame(Gene=character(), animal=factor(), value=numeric())
for (i in 1:nrow(hamp_melted)) {
  if(i %% 2 != 0) {
    j = i+1
    new_val <- mean(c(hamp_melted[i,3],hamp_melted[j,3]))
    hamp_avg <- rbind(hamp_avg, new_val)
  }
}
hamp_avg$pig <- paste("X", c(1:88), sep="")
hamp_avg$gene <- "HAMP"
colnames(hamp_avg) <- c("value", "variable", "Gene")

hamp_merged <- merge(hamp_avg, gapdh, by="variable", suffix=c(".HAMP", ".GAPDH"))
hamp_merged <- hamp_merged[,-4]
hamp_merged$normalized = hamp_merged$value.HAMP - hamp_merged$value.GAPDH
hamp_merged <- hamp_merged[order(hamp_merged["normalized"]),]
hamp_merged$rank <- c(1:nrow(hamp_merged))
hamp_merged1 <- data.frame(Pig = hamp_merged$variable, Gene = hamp_merged$Gene.HAMP, orig_value = hamp_merged$value.HAMP, GAPDH = hamp_merged$value.GAPDH, normalized = hamp_merged$normalized, rank = hamp_merged$rank)


####IMPORT DDX58
ddx <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/ddx58.txt", h=T)
ddx <- melt(ddx)
ddx_avg <- data.frame(Gene=character(), animal=factor(), value=numeric())
for (i in 1:nrow(ddx)) {
  if(i %% 2 != 0) {
    j = i+1
    new_val <- mean(c(ddx[i,3],ddx[j,3]))
    ddx_avg <- rbind(ddx_avg, new_val)
  }
}
ddx_avg$pig <- paste("X", c(1:88), sep="")
ddx_avg$gene <- "DDX58"
colnames(ddx_avg) <- c("value", "variable", "Gene")

ddx_merged <- merge(ddx_avg, gapdh, by="variable", suffix=c(".DDX58", ".GAPDH"))
ddx_merged <- ddx_merged[,-4]
ddx_merged$normalized = ddx_merged$value.DDX58 - ddx_merged$value.GAPDH
ddx_merged <- ddx_merged[order(ddx_merged["normalized"]),]
ddx_merged$rank <- c(1:nrow(ddx_merged))
ddx_merged1 <- data.frame(Pig = ddx_merged$variable, Gene = ddx_merged$Gene.DDX58, orig_value = ddx_merged$value.DDX58, GAPDH = ddx_merged$value.GAPDH, normalized = ddx_merged$normalized, rank = ddx_merged$rank)


####COMBINE INTO A MASTER DATA FRAME
master <- rbind(hamp_merged1, hp_merged, ddx_merged1, nlrp5_merged)

##remove piggies - outliers.
master_trimmed <- subset(master, master$Pig != "X33")
master_trimmed <- subset(master_trimmed, master_trimmed$Pig != "X74")

##ALL PIGS
ggplot(master, aes(x=rank, y=normalized)) + geom_point(shape = 18, size=3) + facet_wrap(~Gene) + theme_bw() + ylab(bquote('Gene expression relative to GAPDH ('~log[2]~')')) + xlab("") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + geom_line() + scale_y_continuous(breaks=seq(-10,10,1)) + geom_hline(yintercept = 0)

##FOR PLOTTING INDIVIDUAL GENES
#ggplot(hamp_sorted, aes(x=rank, y=value)) + geom_point(shape = 18, size=3) + facet_wrap(~Gene) + theme_bw() + ylab(bquote('Gene expression relative to GAPDH ('~log[2]~')')) + xlab("") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + geom_line() + scale_y_continuous(breaks=seq(-10,10,1)) + geom_hline(yintercept = 0)

## ALL OF THE GENES, OUTLIERS REMOVED. 
ggplot(master_trimmed, aes(x=rank, y=normalized)) + 
  geom_point(shape = 18, size=3) + 
  facet_wrap(~Gene) + 
  theme_bw() + 
  ylab(bquote('Gene expression relative to GAPDH ('~log[2]~')')) + 
  xlab("") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line(colour="light grey"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + 
  theme(strip.text = element_text(size=12, face = "italic")) +
  geom_line() + 
  scale_y_continuous(breaks=seq(-10,10,1)) + 
  geom_hline(yintercept=0) #+
  #theme(strip.background = element_rect(colour="black", fill="light grey"))

