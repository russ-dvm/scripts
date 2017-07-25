#RS Fraser
#####NOTE#####
#The first 180 or so lines are the original implementation of the figures. Following discussions with BL on 17-07-25, we determined that a) we should be using the same averages as were used in the heatmap figure; and b) the values should NOT be GAPDH normalized. New figure script starts at the end.

library(reshape2)
library(ggplot2)
library(tidyr)
library(psych)
library(plyr)


###IMPORT REF GENE (GAPDH)
# gapdh <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/gapdh.txt", h=T)
# gapdh <- melt(gapdh)


###IMPORT HAPTOGLOBIN -- multiple probes. Disentangle.
# hp_raw = read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/haptoglobin.txt", h = T)
# hp_raw_melt <- melt(hp_raw)
# hp_raw_melt <- hp_raw_melt[order(hp_raw_melt["value"], decreasing = T), ]
# z <- y
# z$rank <- c(1:nrow(y))

#Check to see if all probes correlate well.
# check <- spread(hp_raw_melt, probe, value)
# corr.test(check[2:6])
###probe 667554 is an outlier.

##subset based on probe
# p1 <- subset(hp_raw_melt, hp_raw_melt$probe == "A_72_P178006")
# p2 <- subset(hp_raw_melt, hp_raw_melt$probe == "A_72_P568194")
# p3 <- subset(hp_raw_melt, hp_raw_melt$probe == "A_72_P581127")
# p4 <- subset(hp_raw_melt, hp_raw_melt$probe == "A_72_P657000")
# p5 <- subset(hp_raw_melt, hp_raw_melt$probe == "A_72_P667554")
# 
# hp_list <- list(p1, p2, p3, p4)

##Normalize the value for each probe
# 
# hp_list_normalized <- lapply(hp_list, merge, gapdh, by.x = "variable", by.y='variable')
# hp_list_normalized <- lapply(hp_list_normalized, function(z) {z$norm <- z$value.x - z$value.y; z})


# complete <- rbind(p1, p2, p3, p4, p5)
# complete <- ldply(hp_list_normalized, rbind)

#average the five probes
# avg <- data.frame(Pig = as.character(), value=as.integer())

#Because all the probes and pigs have their normalized values at this stage, you can choose between averaging the original values, or averaging the gapdh normalized values. "value.x" == the original values; "norm" == gapdh normalized. BL and I did some napkin math and decided that there was no difference between the avg(normalized values) and (avg(original values) - gapdh). 
# for (i in 1:length(unique(complete$variable))){
#   rank <- subset(complete, variable==paste("X", i, sep=""))
#   avg <- rbind(avg, c(rank$variable[1], as.numeric(mean(rank$value.x))))
#   }

##MAIN HP TABLE
# avg$Gene <- "HP"
# colnames(avg) <- c("variable", "value", "Gene")
# hp_sorted <- data.frame(Gene=avg$Gene, variable = avg$variable, value = avg$value)
# hp_sorted$variable = paste("X", hp_sorted$variable, sep = "")
# hp_merged <- merge(hp_sorted, gapdh, by = "variable", suffixes = c(".HP", ".GAPDH"))
# hp_merged <- hp_merged[, -4]
# hp_merged$normalized <- hp_merged$value.HP - hp_merged$value.GAPDH
# 
# hp_merged <- hp_merged[order(hp_merged$normalized),]
# hp_merged$rank <- c(1:nrow(hp_merged))
# colnames(hp_merged) <- c("Pig", "Gene", "orig_value", "GAPDH", "normalized", "rank")

##A lot of time was spent trying to reconcile differences between my graphs and Hein Snyman's original graphs in his thesis. After considerable digging and troubleshooting, looks like the HP graph was not normalized to GAPDH. Here is a version of an un-normalized HP graph for comparison. 
# hein_style_hp <- hp_merged[order(hp_merged$orig_value),]
# hein_style_hp$rank <- c(1:nrow(hein_style_hp))
# ggplot(hein_style_hp, aes(x=rank, y=orig_value)) + geom_point()

####IMPORT MBL2 & NLRP5
# mbl2_raw <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/mbl2.txt", h=T)
# mbl2_melted <- melt(mbl2_raw)
# mbl2_merged <- merge(mbl2_melted, gapdh, by = "variable", suffix = c(".MBL2", ".GAPDH"))
# mbl2_merged <- mbl2_merged[, -4]
# mbl2_merged$normalized <- mbl2_merged$value.MBL2 - mbl2_merged$value.GAPDH
# mbl2_merged <- mbl2_merged[order(mbl2_merged["normalized"]),]
# mbl2_merged$rank <- c(1:nrow(mbl2_merged))
# colnames(mbl2_merged) <- c("Pig", "Gene", "orig_value", "GAPDH", "normalized", "rank")
# 
# 
# nlrp5 <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/nlrp5.txt", h=T)
# nlrp5_melted <- melt(nlrp5)
# nlrp5_merged <- merge(nlrp5_melted, gapdh, by="variable", suffix=c(".NLRP5", ".GAPDH"))
# nlrp5_merged <- nlrp5_merged[,-4]
# nlrp5_merged$normalized <- nlrp5_merged$value.NLRP5 - nlrp5_merged$value.GAPDH
# nlrp5_merged <- nlrp5_merged[order(nlrp5_merged["normalized"]),]
# nlrp5_merged$rank <- c(1:nrow(nlrp5_merged))
# colnames(nlrp5_merged) <- c("Pig", "Gene", "orig_value", "GAPDH", "normalized", "rank")

##HAMP
# hamp_raw <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/hamp.txt", h=T)
# hamp_melted <- melt(hamp_raw)
# hamp_avg <- data.frame(Gene=character(), animal=factor(), value=numeric())
# for (i in 1:nrow(hamp_melted)) {
#   if(i %% 2 != 0) {
#     j = i+1
#     new_val <- mean(c(hamp_melted[i,3],hamp_melted[j,3]))
#     hamp_avg <- rbind(hamp_avg, new_val)
#   }
# }
# hamp_avg$pig <- paste("X", c(1:88), sep="")
# hamp_avg$gene <- "HAMP"
# colnames(hamp_avg) <- c("value", "variable", "Gene")
# 
# hamp_merged <- merge(hamp_avg, gapdh, by="variable", suffix=c(".HAMP", ".GAPDH"))
# hamp_merged <- hamp_merged[,-4]
# hamp_merged$normalized = hamp_merged$value.HAMP - hamp_merged$value.GAPDH
# hamp_merged <- hamp_merged[order(hamp_merged["normalized"]),]
# hamp_merged$rank <- c(1:nrow(hamp_merged))
# hamp_merged1 <- data.frame(Pig = hamp_merged$variable, Gene = hamp_merged$Gene.HAMP, orig_value = hamp_merged$value.HAMP, GAPDH = hamp_merged$value.GAPDH, normalized = hamp_merged$normalized, rank = hamp_merged$rank)


####IMPORT DDX58
# ddx <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/ddx58.txt", h=T)
# ddx <- melt(ddx)
# ddx_avg <- data.frame(Gene=character(), animal=factor(), value=numeric())
# for (i in 1:nrow(ddx)) {
#   if(i %% 2 != 0) {
#     j = i+1
#     new_val <- mean(c(ddx[i,3],ddx[j,3]))
#     ddx_avg <- rbind(ddx_avg, new_val)
#   }
# }
# ddx_avg$pig <- paste("X", c(1:88), sep="")
# ddx_avg$gene <- "DDX58"
# colnames(ddx_avg) <- c("value", "variable", "Gene")
# 
# ddx_merged <- merge(ddx_avg, gapdh, by="variable", suffix=c(".DDX58", ".GAPDH"))
# ddx_merged <- ddx_merged[,-4]
# ddx_merged$normalized = ddx_merged$value.DDX58 - ddx_merged$value.GAPDH
# ddx_merged <- ddx_merged[order(ddx_merged["normalized"]),]
# ddx_merged$rank <- c(1:nrow(ddx_merged))
# ddx_merged1 <- data.frame(Pig = ddx_merged$variable, Gene = ddx_merged$Gene.DDX58, orig_value = ddx_merged$value.DDX58, GAPDH = ddx_merged$value.GAPDH, normalized = ddx_merged$normalized, rank = ddx_merged$rank)

####COMBINE INTO A MASTER DATA FRAME
# master_old <- rbind(hamp_merged1, hp_merged, ddx_merged1, nlrp5_merged)

##remove outlier piggies: X33 = 357, X74 = 812
# master_trimmed_old <- subset(master_old, master_old$Pig != "X33")
# master_trimmed_old <- subset(master_trimmed_old, master_trimmed_old$Pig != "X74")

##ALL PIGS (OUTLIERS PRESENT)
# ggplot(master, aes(x=rank, y=normalized)) + 
#   geom_point(shape = 18, size=3) + 
#   facet_wrap(~Gene) + theme_bw() + 
#   ylab(bquote('Gene expression relative to GAPDH ('~log[2]~')')) + 
#   xlab("") + 
#   theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + 
#   geom_line() + 
#   scale_y_continuous(breaks=seq(-10,10,1)) + 
#   geom_hline(yintercept = 0)

##FOR PLOTTING INDIVIDUAL GENES
# ggplot(hamp_sorted, aes(x=rank, y=value)) + 
#   geom_point(shape = 18, size=3) + 
#   facet_wrap(~Gene) + 
#   theme_bw() + 
#   ylab(bquote('Gene expression relative to GAPDH ('~log[2]~')')) + 
#   xlab("") + 
#   theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text = element_text(face = "italic"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + 
#   geom_line() + 
#   scale_y_continuous(breaks=seq(-10,10,1)) + 
#   geom_hline(yintercept = 0)

### ALL OF THE GENES, OUTLIERS REMOVED. 
# fig2 <- ggplot(master_trimmed, aes(x=rank, y=normalized)) + 
#   geom_line() + 
#   geom_hline(yintercept=0) +
# 
#   geom_point(shape = 23, size=3, colour="black", fill="dark grey") + 
#   facet_wrap(~Gene) + 
#   theme_bw() + 
#   ylab(expression(paste("Gene expression relative to", italic(" GAPDH "), '('~log[2]~')'))) +
#   xlab("") + 
#   theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line(colour="light grey"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + 
#   theme(strip.text = element_text(size=12, face = "italic"), text = element_text(size=12)) +
#   scale_y_continuous(breaks=seq(-10,10,1))
# fig2
# ggsave(fig2, file="~/Desktop/Fig2.eps", width = 174, height = 174, units = "mm", dpi = 600 )

###UPDATED FIGURE SCRIPT 17-07-25
heatmap_data <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/heatmap_data.txt", h=T, sep="\t")

hamp <- heatmap_data[, c(1,grep("HAMP", colnames(heatmap_data)))]
hamp_sorted <- hamp[order(hamp$HAMP.AVG),]
hamp_sorted$rank <- c(1:nrow(hamp_sorted))
hamp_sorted$Gene <- "HAMP"
colnames(hamp_sorted)[2] <- "Value"

hp <- heatmap_data[, c(1,grep("HP.AVG", colnames(heatmap_data)))]
hp_sorted <- hp[order(hp$HP.AVG),]
hp_sorted$rank <- c(1:nrow(hp_sorted))
hp_sorted$Gene <- "HP"
colnames(hp_sorted)[2] <- "Value"

ddx58 <- heatmap_data[, c(1,grep("DDX58", colnames(heatmap_data)))]
ddx58_sorted <- ddx58[order(ddx58$DDX58.AVG),]
ddx58_sorted$rank <- c(1:nrow(ddx58_sorted))
ddx58_sorted$Gene <- "DDX58"
colnames(ddx58_sorted)[2] <- "Value"

nlrp <- heatmap_data[, c(1,grep("NLRP5", colnames(heatmap_data)))]
nlrp_sorted <- nlrp[order(nlrp$NLRP5),]
nlrp_sorted$rank <- c(1:nrow(nlrp_sorted))
nlrp_sorted$Gene <- "NLRP5"
colnames(nlrp_sorted)[2] <- "Value"

#Combine the genes
master <- rbind(hamp_sorted, hp_sorted, ddx58_sorted, nlrp_sorted)
master$Gene <- factor(master$Gene, levels = c("HAMP", "HP", "DDX58", "NLRP5"))

##remove outlier piggies: X33 = 357, X74 = 812
master_trimmed <- subset(master, master$pig != 357)
master_trimmed <- subset(master_trimmed, master_trimmed$pig != 812)

#Plot
fig2 <- ggplot(master_trimmed, aes(x=rank, y=Value)) + 
  geom_line() + 
  geom_hline(yintercept=0) +
  geom_point(shape = 23, size=3, colour="black", fill="dark grey") + 
  facet_wrap(~Gene) + 
  theme_bw() + 
  ylab(expression(paste("Gene expression", '('~log[2]~')'))) +
  xlab("") + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_line(colour="light grey"), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="light grey"), panel.grid.minor.y = element_blank()) + 
  theme(strip.text = element_text(size=12, face = "italic"), text = element_text(size=12)) +
  scale_y_continuous(breaks=seq(-10,10,1))
fig2
ggsave(fig2, file="~/Desktop/fig2.eps", width = 174, height = 174, units = "mm", dpi = 600 )

