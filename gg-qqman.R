library(ggplot2)
library(qqman)
library(viridis)
library(ggforce)


###STEPHEN TURNER'S CODE FROM 
## https://github.com/stephenturner/qqman/blob/master/R/manhattan.R

manhattan <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", 
                      col=c("gray10", "gray60"), chrlabs=NULL,
                      suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), 
                      highlight=NULL, logp=TRUE, annotatePval = NULL, annotateTop = TRUE, ...) {
  
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
    xlabel = paste('Chromosome',unique(d$CHR),'position')
    labs = ticks
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
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
      
    }
    xlabel = 'Chromosome'
    #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
    labs <- unique(d$CHR)
  }
  
  # Initialize plot
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  
  #Direct ticks and labels into their own values and return the data frame - ggplot addition at end (outside the function)
  assign('ticks', ticks, envir=.GlobalEnv)
  assign('labs',labs, envir= .GlobalEnv)
  return(d)
  
}

####IMPORT DATA

##EQUINE
results <- read.table('~/playground/qqman/results2.txt', h=T)

##BOVINE
#results <- read.table("", h=T)


##Popoolation2 P values have already been -log, qqman expects raw P values, so have to re-calculate raw P
results$P <- 10^(-results$P)

### Use qqman manhattan function to format the data appropriately for a Manhattan plot
d <- manhattan(results)

###Determine BH significance levels###
#order data by p value, from lowest to highest
sorted.by.p <- d[order(d$P),]	
#determine number of rows
rownum <- nrow(sorted.by.p) 											
#assign ranks to each row, lowest p value = 1
sorted.by.p$rank <- c(1:rownum) 										
#set the false discovery rate
fdr = 0.01																
#calculate the benjamini hochberg critical value - the integer = False Discovery Rate
benj.hoch <- (sorted.by.p$rank/rownum)*fdr 		
#total #in n_D.fet = 1075951

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

#inverse log of the p-value of highest ranked sig value
intercept <- lowest.significant$P				

##Change the chroms from numeric to factors, for colouring purposes.
d$CHR <- as.factor(d$CHR)


###"extra" data points that were added to set the pos values should be removed prior to plotting:
#Requires the snp value "hide" for any additional "data" point that was added.
d.real <- subset(d, d$SNP != 'hide')


##GGPLOT MOD BY RUSSELL FRASER
ggplot(d.real) + geom_point(aes(x=pos, y=logp, color=CHR)) + theme_classic()+ ylab("-log(p)") + theme(axis.text.x=element_text(angle = 60, hjust=1, size=9), legend.title=element_blank(), text=element_text(size=20), plot.title=element_text(size=15), legend.position = "none") + scale_color_viridis(discrete=T) + geom_hline(yintercept=-log10(intercept), colour="red") + facet_zoom(x=CHR==1)+ facet_zoom(x=CHR==1 & BP > 88892979 & BP < 89006643)# + scale_x_continuous(breaks=c(ticks), labels=labs)

# x=CHR==1 & BP > 88892979 & BP < 89006643
# x=CHR==8 & BP > 40893025 & BP < 41121736

#scale_color_manual(values=c("orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange", "navy blue","orange"))

