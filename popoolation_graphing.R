library(ggplot2)
library(viridis)
library(ggforce)
#BOVINE
#popoolation <- read.table("~/bovine/merged_runs/popoolation/subsampled200-no-na.fet", col.names=c("Chrom", "Position", "A", "B", "C", "Pop", "invlogp"))

#EQUINE
popoolation <- read.table("~/equine/2014_11_24/popoolation/subsampled200-tabbed.fet", col.names=c("Chrom", "Position", "A", "B", "C", "Pop", "invlogp"))


###Determine BH significance levels

#order data by p value, from lowest to highest
sorted.by.p <- popoolation[order(popoolation$invlogp, decreasing=T),]	

#determine number of rows
rownum <- nrow(sorted.by.p) 											

#assign ranks to each row, lowest p value = 1
sorted.by.p$rank <- c(1:rownum) 										

#reverse the invlogp calculation
sorted.by.p$actual.p <- 10^(-sorted.by.p$invlogp)

#set the false discovery rate
fdr = 0.01																

#calculate the benjamini hochberg critical value - the integer = False Discovery Rate
benj.hoch <- (sorted.by.p$rank/rownum)*fdr 		
#total #in n_D.fet = 1075951

#add benjamini hochbergs to the data frame
sorted.by.p$bh <- benj.hoch 							

#determine significant (i.e. is the p value less than the BH critical value?)
sig <- sorted.by.p$actual.p < sorted.by.p$bh 		

#add the significance to the data frame
sorted.by.p$significant <- sig 				

#convert the continuous integer CHR to a factor for colour.
Chromosome <- as.factor(sorted.by.p$Chrom)		

#add it to the data frame
sorted.by.p$Chromosome <- Chromosome									

###Calculate the intercept line (-logp of the highest p < bh)

#Take a subset composed only of values that are significant
all.significant <- subset(sorted.by.p, significant == TRUE)		

#Take the last significant one (i.e. highest rank)
lowest.significant <- tail(all.significant, 1)	
lowest.significant
lowest.significant.row <- lowest.significant$rank

#inverse log of the p-value of highest ranked sig value
intercept <- lowest.significant$invlogp								

#BOVINE - order by chromosome number
#sorted.by.p$Chromosome <- factor(sorted.by.p$Chromosome, levels = c("chr1", "chr2", "chr3", "chr7", "chr8", "chr11", "chr14", "chr16", "chr18", "chr24", "chr26", "chr28"))

#EQUINE - order by chromosome number
sorted.by.p$Chromosome <- factor(sorted.by.p$Chromosome, levels = c("chr1", "chr2", "chr5", "chr8", "chr9", "chr10", "chr15", "chr19", "chr21", "chr25"))



plot1 <- ggplot(sorted.by.p) + geom_jitter(aes(x=Chromosome, y=invlogp, color=Chromosome), width=0.5) + geom_hline(yintercept=intercept, col="red") + theme_classic() +  xlab("")+ ylab("-log(p)")+ theme(legend.position="none") + expand_limits(y=c(0,6)) + scale_colour_manual(values=c("orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue")) + theme(axis.text.x=element_text(angle = 45, hjust=1), legend.title=element_blank(), text=element_text(size=20), plot.title=element_text(size=15))
plot1

sig_only <- subset(sorted.by.p, sorted.by.p$significant == "TRUE")
sig_only <- head(sorted.by.p, paste(lowest.significant.row))
nrow(sig_only)
write.table(sig_only, file = "~/Desktop/sig-only200.txt", quote = F, sep = "\t", row.names = F)


  # ggd.qqplot = function(pvector, main=NULL, ...) {
#   o = -log10(sort(pvector,decreasing=F))
#   e = -log10( 1:length(o)/length(o) )
#   plot(e,o,pch=19,cex=1, main=main, ...,
#        xlab=expression(Expected~~-log[10](italic(p))),
#        ylab=expression(Observed~~-log[10](italic(p))),
#        xlim=c(0,max(e)), ylim=c(0,max(o)))
#   lines(e,e,col="red")
# }
# 
# # Using the ggd.qqplot() function
# ggd.qqplot(sorted.by.p$actual.p)
# 
# # Add a title
# ggd.qqplot(sorted.by.p$actual.p, "QQ-plot of p-values using ggd.qqplot")
# nrow(sorted.by.p)
# 
# hist(sorted.by.p$invlogp)
# hist(sorted.by.p$actual.p)

#######################
###ZOOM IN ON REGION###
#######################

ggplot(sorted.by.p) + geom_jitter(aes(x=Chromosome, y=invlogp, color=Chromosome), width=0.49) + facet_zoom(x=Chromosome=='chr10') + theme_bw() + xlab("") + ylab('Inv log (P)') + theme(legend.position="none") + scale_color_viridis(discrete=T)


+ geom_hline(yintercept=intercept, col="red") + theme_classic() +  xlab("")+ ylab("-log(p)")+ theme(legend.position="none") + expand_limits(y=c(0,6)) + 
  