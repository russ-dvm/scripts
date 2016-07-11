#Calculate the fisher exact test value on the output of 
#fisher-test.pl  Input might be a vcf file.
#RSF, rfrase03@uoguelph.ca
#January 11, 2016


#packages
library(ggplot2)

#input data
fet <- read.table("~/sandbox/allele_freq_comparison/gatk_depths/test.vcf", h=T)

#set empty variable containers
nref=NULL
nalt=NULL
dref=NULL
dalt=NULL
pvalue=NULL
logp=NULL
ref.freq=NULL
alt.freq=NULL

#calculate fishers exact test statistic and the inverse log and the allele frequencies for ref and alt
for (i in 1:nrow(fet)) {
	nref[i] <- fet[i,8]
	nalt[i] <- fet[i,9]
	dref[i] <- fet[i,10]
	dalt[i] <- fet[i,11]
	pvalue[i] <- fisher.test(rbind(c(nref[i],dref[i]), c(nalt[i],dalt[i])))$p.value
	logp[i] <- -log10(pvalue[i])
	ref.freq[i] <- nalt[i]/(nref[i]+nalt[i])
	alt.freq[i] <- dalt[i]/(dref[i]+dalt[i])
	}

#append the data to the original table	
fet$pvalue <- pvalue
fet$logp <- logp
fet$ref.freq <- ref.freq
fet$alt.freq <- alt.freq


#do some benjamini-hochberg stuff
sorted.by.p <- fet[order(fet$pvalue, decreasing=F),]						#order data by p value, from lowest to highest
rownum <- nrow(sorted.by.p) 											#determine number of rows
sorted.by.p$rank <- c(1:rownum) 										#assign ranks to each row, lowest p value = 1
fdr = 0.01																#set the false discovery rate
benj.hoch <- (sorted.by.p$rank/rownum)*fdr 								#calculate the benjamini hochberg critical value - the integer = False Discovery Rate
sorted.by.p$bh <- benj.hoch 											#add benjamini hochbergs to the data frame
sig <- sorted.by.p$pvalue < sorted.by.p$bh 									#determine significant (i.e. is the p value less than the BH critical value?)
sorted.by.p$significant <- sig 											#add the significance to the data frame
Chromosome <- as.factor(sorted.by.p$chrom)								#convert the continuous integer CHR to a factor for colour.
sorted.by.p$Chromosome <- Chromosome									#add it to the data frame

#Calculate the intercept line (-logp of the highest p < bh)
all.significant <- subset(sorted.by.p, significant == TRUE)				#Take a subset composed only of values that are significant
lowest.significant <- tail(all.significant, 1)							#Take the last significant one (i.e. highest rank)
intercept <- lowest.significant$logp								#inverse log of the p-value of highest ranked sig value


##SUBSETS##
passed <- subset(sorted.by.p, filter=="PASS")
chr1 <- subset(sorted.by.p, chrom=="chr1")
chr1.spd <- subset(chr1, position > 8e7)




###OUTPUT###
plot1 <- ggplot(sorted.by.p) + 
	geom_point(aes(x=Chromosome, y=logp, colour=Chromosome), position=position_jitter(width=1)) + 			#modify the width values to separate points for each chr
	geom_hline(yintercept=intercept, col="red") + theme_classic() + 															#red significance line - can probably change to dashed, etc
	xlab("")+																							#x-axis label
	ylab("-log(p)")+
	theme(legend.position="none")  																			#suppress the legend



plot1

#output the data
write.table(sorted.by.p, file="~/sandbox/allele_freq_comparison/gatk_depths/fisher-exact-test-stats.txt", row.names=F, sep="\t", quote=F)
ggsave(plot = plot1, filename="~/fisher_output/plot1.png")