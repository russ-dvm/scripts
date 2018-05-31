library(tidyverse)
library(qqman)
library(viridis)
library(ggforce)


######ASSESS THE DATA
##IMPORT

freq_table <- read.table("~/equine/2014_11_24/popoolation/subsampled370.17-08-31.freqs.no-plgyrps_rc", h=T, stringsAsFactors = F)
fisher <- read.table("~/equine/2014_11_24/popoolation/subsampled370.17-08-31.ready.no-pglyrp.fet", col.names=c("Chrom", "pos", "A", "B", "C", "Pop", "invlogp"))

##Check the the major and minor alleles for each population are the same
freq_table <- separate(freq_table, major_alleles.maa., c("maa_n", "maa_d"), sep = 1)
freq_table <- separate(freq_table, minor_alleles.mia., c("mia_n", "mia_d"), sep = 1)
freq_table$maa_same <- freq_table$maa_n == freq_table$maa_d
freq_table$mia_same <- freq_table$mia_n == freq_table$mia_d


##Fix table where major/minor alleles are different for the two pops
for (i in 1:nrow(freq_table)){
  if (freq_table[i,]$maa_same == F){
    storage_maa_d <- freq_table[i,]$maa_d
    storage_maa_2 <- freq_table[i,]$maa_2
    freq_table[i,]$maa_d <- freq_table[i,]$mia_d
    freq_table[i,]$maa_2 <- freq_table[i,]$mia_2
    freq_table[i,]$mia_d <- storage_maa_d
    freq_table[i,]$mia_2 <- storage_maa_2
  }
  
}
##Re-check that the values are the same
freq_table$maa_same <- freq_table$maa_n == freq_table$maa_d
freq_table$mia_same <- freq_table$mia_n == freq_table$mia_d

##Good idea to inspect the table manually -- for the equine data there are a few SNPs that should be excluded and a couple that need to be manually fixed.

##Have a look at multi-allelic sites:
multi_allele <- subset(freq_table, nchar(freq_table$allele_states) > 3)
##How do those look in the FET data and raw data? 
## Does it matter? These alleles were excluded in SNV calling; should be excluded here too.

##Fix things that are bad
freq_table[1763,]$mia_d <- "A"
freq_table[19,]$mia_d <- "A"
freq_table$maa_same <- freq_table$maa_n == freq_table$maa_d
freq_table$mia_same <- freq_table$mia_n == freq_table$mia_d

##Subset with only passing comparisons.
freq_table <- freq_table[freq_table$mia_same == T,]

##Merge the results of fisher's exact test with the actual frequencies
mega_table <- merge(freq_table, fisher, by.x = "pos")
mega_table_curated <- mega_table[c(-3, -4, -6, -7, -16:-22)]


####Manipulate the format of frequency diffs
##Isolate the minor allele fractions
maybe <- data.frame("mia_n" = as.character(mega_table$mia_1), "mia_d" = as.character(mega_table$mia_2), "pos" = as.numeric(mega_table$pos), "chrom" = as.character(mega_table$chr))
##Split numerator and denominator
maybe <- separate(maybe, mia_n, c("alt_n", "total_n"), sep = "/")
maybe <- separate(maybe, mia_d, c("alt_d", "total_d"), sep = "/")
##Convert back to numbers
maybe$alt_n <- as.numeric(maybe$alt_n)
maybe$alt_d <- as.numeric(maybe$alt_d)
maybe$total_n <- as.numeric(maybe$total_n)
maybe$total_d <- as.numeric(maybe$total_d)
##Calculate the frequency
maybe$n_freq <- maybe$alt_n/maybe$total_n
maybe$d_freq <- maybe$alt_d/maybe$total_d
##Calculate the difference in the frequencies -- SHOULD THIS A DIVISION?? 
maybe$diff <- maybe$n_freq - maybe$d_freq

##Do some QC - have a look at the distribution of differences
library(MASS)
ggplot(maybe, aes(x = diff)) + geom_histogram(colour = "black", fill = "white") + theme_bw()
qqnorm(maybe$diff) + qqline(maybe$diff)


###Merge back into main table and create a final publication ready table
mega_table_curated <- merge(mega_table_curated, maybe, by.x = "pos")
final_table <- mega_table_curated[,c(2,1,4,6,10,11,18:20,12)]

colnames(final_table) <- c("Chrom", "Position", "Ref allele", "Alt allele", "Frac (reads) in normal", "Frac (reads) in diseasd", "Freq in normal", "Freq in diseased", "Freq diff", "-log(p)")

##Sort by p value and assign rank
final_table <- final_table[order(final_table$`-log(p)`, decreasing = T),]
final_table$rank <- c(1:nrow(final_table)) 

##Print final table
write.table(final_table, "~/equine/2014_11_24/popoolation/pub_table.txt", sep="\t", row.names=F, quote=F)
write.table(final_table, "~/Dropbox/temp/pub_table.txt", sep="\t", row.names=F, quote=F)

########GRAPHING########

###STEPHEN TURNER'S CODE FROM 
## https://github.com/stephenturner/qqman/blob/master/R/manhattan.R
manhattan <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", 
                      col=c("gray10", "gray60"), chrlabs=NULL,
                      suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), 
                      highlight=NULL, logp=TRUE, annotatePval = NULL, annotateTop = TRUE, X=NULL, Y=NULL, ...) {
  
  # Not sure why, but package check will warn without this.
  CHR=BP=P=index=NULL
  
  # Check for sensible dataset
  ## Make sure you have chr, bp and p columns.
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
  ## warn if you don't have a snp column
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  ## make sure chr, bp, and p columns are numeric.
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  
  # Create a new data.frame with columns called CHR, BP, and P.
  d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]])
  
  # If the input data frame has a SNP column, add it to the new data frame you're creating.
  if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
  
  # Set positions, ticks, and labels for plotting
  ## Sort and keep only values where is numeric.
  #d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  #d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
  if (logp) {
    d$logp <- -log10(d$P)
  } else {
    d$logp <- d$P
  }
  d$pos=NA
  
  
  # Fixes the bug where one chromosome is missing by adding a sequential index column.
  d$index=NA
  ind = 0
  for (i in unique(d$CHR)){
    ind = ind + 1
    d[d$CHR==i,]$index = ind
  }
  
  # This section sets up positions and ticks. Ticks should be placed in the
  # middle of a chromosome. The a new pos column is added that keeps a running
  # sum of the positions of each successive chromsome. For example:
  # chr bp pos
  # 1   1  1
  # 1   2  2
  # 2   1  3
  # 2   2  4
  # 3   1  5
  nchr = length(unique(d$CHR))
  if (nchr==1) { ## For a single chromosome
    ## Uncomment the next two linex to plot single chr results in Mb
    #options(scipen=999)
    #d$pos=d$BP/1e6
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
    minor_ticks = floor(length(d$pos))/2+1
    xlabel = paste('Chromosome',unique(d$CHR),'position')
    labs = ticks
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
    minor_ticks=NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
        d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
      }
      # Old way: assumes SNPs evenly distributed
      # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
      # New way: doesn't make that assumption
      ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
      
      ##MOD: ticks that bound the chromosome, rather than appear at its center - RSF
      minor_ticks = c(minor_ticks, min(d[d$index == i,]$pos, max(d[d$index == i,]$pos)))
      
    }
    xlabel = 'Chromosome'
    #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
    labs <- unique(d$CHR)
    
    ##MOD - use the X and Y options to rename the X and Y chromosomes to the more appropriate "X" and "Y". This does coerce labs into character instead of numeric but that doesn't seem to affect the output of ggplot2. RSF.
    labs[X] <- "X"
    labs[Y] <- "Y"
  }
  
  # Initialize plot
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  
  ###MODIFIED HERE -- Direct ticks and labels into their own values and return the data frame 
  ## The datafrmae (d) can then be used in ggplot along with the ticks and labels
  assign('ticks', ticks, envir=.GlobalEnv)
  assign('labs', labs, envir= .GlobalEnv)
  assign('minor_ticks', minor_ticks, envir= .GlobalEnv)
  return(d)
  
}

###Import start/stop chromosome positions for equine
start_stop <- read.table("~/equine/2014_11_24/popoolation/start_stop_for_qqman.txt", sep = "\t", h = T)

## Reformat for qqman -- expects three columns, Chr, Bp, P, and SNP-id
man_frame <- data.frame("CHR" = final_table$Chrom, "BP" = as.integer(final_table$Position), "P" = 10^(-final_table$`-log(p)`), "SNP" = rep("snp", nrow(final_table)))
man_frame$CHR <- gsub("chr", "", man_frame$CHR)
man_frame$CHR <- as.numeric(man_frame$CHR)

##Add in start and stop "snps" - pvalue 0 - and snp value of hide. This is necessary so that the chr widths are accurate
man_frame <- rbind(man_frame, start_stop)

### Use qqman manhattan function to format the data appropriately for a Manhattan plot
d <- manhattan(man_frame)

###Determine BH significance levels###
#order data by p value, from lowest to highest
sorted.by.p <- d[order(d$P),]	
#determine number of rows
rownum <- nrow(sorted.by.p) 											
#assign ranks to each row, lowest p value = 1
sorted.by.p$rank <- c(1:rownum) 										
#set the false discovery rate
fdr = 0.000001																
#calculate the benjamini hochberg critical value
benj.hoch <- (sorted.by.p$rank/rownum)*fdr 		
#add benjamini hochbergs to the data frame
sorted.by.p$bh <- benj.hoch 							
#determine significant (i.e. is the p value less than the BH critical value?)
sig <- sorted.by.p$P < sorted.by.p$bh 		
#add the significance to the data frame
sorted.by.p$significant <- sig 				

###Calculate the intercept line (-logp of the highest p < bh)
#Take a subset composed only of values that are significant
all.significant <- subset(sorted.by.p, significant == TRUE)		

#Take the last significant one (i.e. highest rank)
lowest.significant <- tail(all.significant, 1)	
lowest.significant.row <- lowest.significant$rank

#Make sure all rows above the lowest row are considered significant
sorted.by.p[1:lowest.significant.row,]$significant <- TRUE

#inverse log of the p-value of highest ranked sig value
intercept <- lowest.significant$P				


## Add the R padjust to the data frame
d$rp <- p.adjust(d$P, method = "BH")

####Change the chroms from numeric to factors, for colouring purposes.
d$CHR <- as.factor(d$CHR)

## Collectin locus gene coords
colec <- c(88892979, 88922677,88925052,88932901,88935579,88947828,88956555, 89006643)

##GGPLOT MOD BY RUSSELL FRASER
#scale_x_continous seems to break the x-axis of the facet_zoom
#There's probably a smart workaround, but for now, the best I can come up with is saving a copy of the graph with
# and without the sacle_x_continous as SVG and then merging the two. 

fig <- ggplot(d) + 
  geom_point(aes(x=pos, y=-log10(rp), color=CHR)) + 
  theme_bw() +
  theme(panel.grid.minor.x = element_line(colour = "light grey"), panel.grid.major.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  ylab("-log(p)") + 
  xlab("") +
  theme(axis.text.x=element_text(angle = 60, hjust=1), legend.title=element_blank(), legend.position = "none") + 
  theme(text = element_text(size = 10)) +
  scale_colour_manual(values = eg_col) + 
  geom_hline(yintercept=-log10(fdr), colour="red") + 
  facet_zoom(x=CHR==1)+ facet_zoom(x=CHR==1 & BP > 88892979 & BP < 89006643)

fig_a <- fig + 
  scale_x_continuous(minor_breaks = minor_ticks, breaks = ticks, labels = labs)
  
fig_b <- fig +
  scale_x_continuous(minor_breaks = colec, breaks = colec, labels = c("88892979", "SFTPA1", "", "MBL1", "", "SFTPD", "", "89006643"))

fig_a
fig_b

ggsave("~/Dropbox/temp/figure_4a.eps", fig_a, dpi = 1000, units = "mm", width = 190, height = 100)
ggsave("~/Dropbox/temp/figure_4b.eps", fig_b, dpi = 1000, units = "mm", width = 190, height = 100)


sig_only <- subset(sorted.by.p, sorted.by.p$significant == T)
nrow(sig_only)

# ggsave("~/Dropbox/chrom.svg")
# scale_x_continuous(breaks=c(ticks), labels=labs) +
write.table(d, "~/Dropbox/temp/tmp.txt", row.names = F, sep = "\t", quote = F)







eg_col <- c(
  "grey27", 
  "grey67",
  "grey27", 
  "grey67",
  "grey27", 
  "grey67",
  "grey27", 
  "grey27",
  "grey67", 
  "grey67",
  "grey27", 
  "grey67",
  "grey27", 
  "grey67",
  "grey27", 
  "grey67",
  "grey27", 
  "grey67",
  "grey67", 
  "grey67",
  "grey27", 
  "grey67",
  "grey27", 
  "grey67",
  "grey27", 
  "grey27",
  "grey27", 
  "grey67",
  "grey27", 
  "grey67",
  "grey27")
