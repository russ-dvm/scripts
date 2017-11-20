library(tidyverse)
library(qqman)
library(viridis)
library(ggforce)


######ASSESS THE DATA
##IMPORT

freq_table <- read.table("~/bovine/merged_runs/popoolation/nd.freqs_rc", h=T, stringsAsFactors = F)
fisher <- read.table("~/bovine/merged_runs/popoolation/nd.subsampled400.ready.fet", col.names=c("Chrom", "pos", "A", "B", "C", "Pop", "invlogp"))

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
freq_table <- subset(freq_table, nchar(freq_table$allele_count) <= 3)

##Subset with only passing comparisons - there are a number that have Ns instead of an NT.
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

ggplot(maybe, aes(x = diff)) + geom_histogram(colour = "black", fill = "white", binwidth = 0.01) + theme_bw()
qqnorm(maybe$diff) + qqline(maybe$diff)


###Merge back into main table and create a final publication ready table
mega_table_curated <- merge(mega_table_curated, maybe, by.x = "pos")
final_table <- mega_table_curated[,c(2,1,4,6,10,11,18:20,12)]

colnames(final_table) <- c("Chrom", "Position", "Ref allele", "Alt allele", "Frac (reads) in normal", "Frac (reads) in diseasd", "Freq in normal", "Freq in diseased", "Freq diff", "-log(p)")

##
final_table$P <- 10^-final_table$`-log(p)`
final_table$adj <- p.adjust(final_table$P, method = "BH")
final_table$logBH <- -log10(final_table$adj)

##Sort by p value and assign rank
final_table <- final_table[order(final_table$logBH, decreasing = T),]
final_table$rank <- c(1:nrow(final_table)) 

### Add in rsIDs
vcf <- read.table("~/bovine/merged_runs/6_variants/lectin/17-11-09.all.lectin.vcf", h = F)
vcf <- vcf[, c(1:3)]
colnames(vcf) <- c("chr", "Position", "rsID")

final_table <- merge(final_table, vcf, by.x = "Position", all.x = T)

##Print final table
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
start_stop <- read.table("~/bovine/merged_runs/popoolation/start_stop_for_qqman.txt", sep = "\t", h = T)

## Reformat for qqman -- expects three columns, Chr, Bp, P, and SNP-id
man_frame <- data.frame("CHR" = final_table$Chrom, "BP" = as.integer(final_table$Position), "P" = 10^(-final_table$`-log(p)`), "SNP" = rep("snp", nrow(final_table)))
man_frame$CHR <- gsub("chr", "", man_frame$CHR)
man_frame$CHR <- as.numeric(man_frame$CHR)

##Add in start and stop "snps" - pvalue 0 - and snp value of hide. This is necessary so that the chr widths are accurate
man_frame <- rbind(man_frame, start_stop)

### Use qqman manhattan function to format the data appropriately for a Manhattan plot
d <- manhattan(man_frame)

###Determine BH significance levels###

d$bh <- p.adjust(d$P, method = "BH")
d$logbh <- -log10(d$bh)
class1 <- 2e-7
class2 <- 5e-6
class3 <- 1e-4

####Change the chroms from numeric to factors, for colouring purposes.
d$CHR <- as.factor(d$CHR)

##Change chr30 to chrX
labs[30] <- "X"
## Collectin locus gene coords
colec <- c(35541900,35591906,35602846,35668986,35700161,35714667,35748478,35813760,35824383,35840849,35846070,35850665,35854765,35870565)
## Adjust so that it reflects the pos argument from the manhattan function. Can determine the distance by pos - BP - then add that value to the above.
colec <- 2414264736 + colec


## No of sig snps
sigAtClass1 <- subset(d, bh <= class1)
sigAtClass2 <- subset(d, bh <= class2)
sigAtClass3 <- subset(d, bh <= class3)

##GGPLOT MOD BY RUSSELL FRASER
#scale_x_continous seems to break the x-axis of the facet_zoom
#There's probably a smart workaround, but for now, the best I can come up with is saving a copy of the graph with
# and without the sacle_x_continous as SVG and then merging the two. 

fig <- ggplot(d) + 
  geom_point(aes(x=pos, y=logbh, color=CHR)) + 
  theme_bw() +
  theme(panel.grid.minor.x = element_line(colour = "light grey"), panel.grid.major.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  # ylab("-log(adjusted italic(p)-value)") +
  ylab(expression(paste("-log(adjusted ", italic("p"), "-value)")))+
  xlab("") +
  theme(axis.text.x=element_text(angle = 60, hjust=1), legend.title=element_blank(), legend.position = "none") + 
  theme(text = element_text(size = 10)) +
  # scale_colour_manual(values = g_col) +
  geom_hline(yintercept=-log10(class1), colour="dark grey", linetype = "dashed") + 
  geom_hline(yintercept=-log10(class2), colour="dark grey") +
  # geom_hline(yintercept=-log10(class3), colour = "dark grey") +
  facet_zoom(x=CHR==28)+ facet_zoom(x=CHR==28 & BP > 35541900 & BP < 35870565) +
  annotate("text", x = 2650906405, y = -log10(class2), label = paste("n =", nrow(sigAtClass2)), vjust = -0.2) +
  annotate("text", x = 2650906405, y = -log10(class1), label = paste("n =", nrow(sigAtClass1)), vjust = -0.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 9))
fig

fig_a <- fig + 
  scale_x_continuous(minor_breaks = minor_ticks, breaks = ticks, labels = labs)

fig_b <- fig +
  scale_x_continuous(minor_breaks = colec, breaks = colec, labels = c(35541900, "CGN", "", "CL46", "", "CL43", "", "SFTPD", "", "MBL1", "", "SFTPA1", "",35870565))

fig_a
fig_b

ggsave("~/Dropbox/temp/figure_4a.eps", fig_a, dpi = 1000, units = "mm", width = 190, height = 100)
ggsave("~/Dropbox/temp/figure_4b.eps", fig_b, dpi = 1000, units = "mm", width = 190, height = 100)



g_col <- c(
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
