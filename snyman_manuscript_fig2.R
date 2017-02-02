library(reshape)
library(ggplot2)
x = read.table("~/Dropbox/temp/haptoglobin.txt", h = T)

y <- melt(x)

y <- y[order(y["value"], decreasing = T), ]
z <- y
z$rank <- c(1:nrow(y))
ggplot(z, aes(x = variable, y = value)) + geom_boxplot() + theme_classic() + xlab("") 
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

#all data
ggplot(complete) + geom_point(aes(x=rank, y=value)) + facet_grid(.~probe) + ggtitle("all probes, all pigs")

#drop the lowest two
ggplot(subset(complete, rank != 1 & rank != 2)) + geom_point(aes(x=rank, y=value)) + facet_grid(.~probe) + ggtitle("all probes, bottom 2 pigs removed")

#average the five probes
avg <- data.frame(Pig = as.character(), value=as.integer())
total = 0
for (i in 1:nrow(p1)){
  rank <- subset(complete, rank==i)
  avg <- rbind(avg, c(rank$variable[1], as.numeric(mean(rank$value))))
  }

avg <- avg[order(avg[,2]),]
avg$rank <- c(1:nrow(avg))

ggplot(avg, aes(x=rank, y=avg[,2])) + geom_point() + ggtitle("average of all 5 probes, all pigs")

avg_no_out <- subset(avg, rank != 1 & rank != 2)
ggplot(avg_no_out) + geom_point(aes(x=rank, y=avg_no_out[,2]), shape=5)+ ggtitle("average of all 5 probes, bottom 2 pigs removed")

#2x2
pub <- rbind(p1, p2, p3, p4)
ggplot(pub, aes(x=rank, y=value)) + geom_point() + facet_wrap(~probe, scales="free") + theme_bw() + ylab("GAPDH") + xlab("Pig ID")

##
avg$Gene <- "HP"
colnames(avg) <- c("variable", "value", "rank", "Gene")
hp_sorted <- data.frame(avg$Gene, avg$variable, avg$value, avg$rank)
colnames(hp_sorted) <- c("Gene", "variable", "value", "rank")
hp_sorted$variable <- paste("X", hp_sorted$variable, sep="")


####IMPORT MBL2 & NLRP5

mbl2_raw <- read.table("~/Dropbox/temp/mbl2.txt", h=T)
mbl2_melted <- melt(mbl2_raw)
mbl2_sorted <- mbl2_melted[order(mbl2_melted["value"]),]
mbl2_sorted$rank <- c(1:nrow(mbl2_sorted))

nlrp5 <- read.table("~/Dropbox/temp/nlrp5.txt", h=T)
nlrp5_melted <- melt(nlrp5)
nlrp5_sorted <- nlrp5_melted[order(nlrp5_melted["value"]),]
nlrp5_sorted$rank <- c(1:nrow(nlrp5_sorted))


##HAMP

hamp_raw <- read.table("~/Dropbox/temp/hamp.txt", h=T)
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
colnames(hamp_avg) <- c("value", "pig", "Gene")
hamp_sorted <- hamp_avg[order(hamp_avg["value"]),]
hamp_sorted$rank <- c(1:nrow(hamp_sorted))
hamp_sorted <- data.frame(hamp_sorted$Gene, hamp_sorted$pig, hamp_sorted$value, hamp_sorted$rank)
colnames(hamp_sorted) <- c("Gene", "variable", "value", "rank")


master <- rbind(hamp_sorted, hp_sorted, mbl2_sorted, nlrp5_sorted)

##remove piggies
master_trimmed <- subset(master, master$variable != "X33")
master_trimmed <- subset(master_trimmed, master_trimmed$variable != "X74")

##ALL PIGS
ggplot(master, aes(x=rank, y=value)) + geom_point(shape = 18, size=3) + facet_wrap(~Gene) + theme_bw() + ylab(bquote('Gene expression relative to GAPDH ('~log[2]~')')) + xlab("") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + geom_line() + scale_y_continuous(breaks=seq(-10,10,1)) + geom_hline(y=0)

##FOR PLOTTING INDIVIDUAL GENES
ggplot(hamp_sorted, aes(x=rank, y=value)) + geom_point(shape = 18, size=3) + facet_wrap(~Gene) + theme_bw() + ylab(bquote('Gene expression relative to GAPDH ('~log[2]~')')) + xlab("") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + geom_line() + scale_y_continuous(breaks=seq(-10,10,1)) + geom_hline(y=0)

## ALL OF THE GENES. 
ggplot(master_trimmed, aes(x=rank, y=value)) + geom_point(shape = 18, size=3) + facet_wrap(~Gene) + theme_bw() + ylab(bquote('Gene expression relative to GAPDH ('~log[2]~')')) + xlab("") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line(colour="light grey"), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + geom_line() + scale_y_continuous(breaks=seq(-10,10,1)) + geom_hline(y=0)

