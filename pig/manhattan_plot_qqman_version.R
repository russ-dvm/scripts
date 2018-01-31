library(ggplot2)
library(qqman)
library(viridis)
library(ggforce)


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

####IMPORT DATA

##eQTL RESULTS
complete_results <-read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/eQTL_results_no-outliers/eQTL_results_LR_min3/cis_summary_86000.txt", h=T)

complete_results$snp_chrom <- sub("chr", "", complete_results$snp_chrom)
complete_results$snp_chrom <- sub("X", "19", complete_results$snp_chrom)

#Determine significance cutoff
sig <- complete_results$p.value < complete_results$FDR



qq_results <- data.frame("CHR" = as.integer(complete_results$snp_chrom), "BP" = as.integer(complete_results$snp_pos), "P" = as.numeric(complete_results$FDR), "SNP" = as.character(complete_results$snp_id))

##In order for this hack to work properly, "data" points corresponding to the first and last position of each chromosome need to be added in. I add them with the P-value as NA, which automatically omits them in ggplot2. Data is from Ensembl Sus scrofa 10.2.
for (x in 1:20){
  qq_results <- rbind(qq_results, c(x, 1, NA, NA))
}
qq_results <- rbind(qq_results, c(1, 315321322, NA, NA))
qq_results <- rbind(qq_results, c(2, 162569375, NA, NA))
qq_results <- rbind(qq_results, c(3, 144787322, NA, NA))
qq_results <- rbind(qq_results, c(4, 143465943, NA, NA))
qq_results <- rbind(qq_results, c(5, 111506441, NA, NA))
qq_results <- rbind(qq_results, c(6, 157765593, NA, NA))
qq_results <- rbind(qq_results, c(7, 134764511, NA, NA))
qq_results <- rbind(qq_results, c(8, 148491826, NA, NA))
qq_results <- rbind(qq_results, c(9, 153670197, NA, NA))
qq_results <- rbind(qq_results, c(10, 79102373, NA, NA))
qq_results <- rbind(qq_results, c(11, 87690581, NA, NA))
qq_results <- rbind(qq_results, c(12, 63558571, NA, NA))
qq_results <- rbind(qq_results, c(13, 218635234, NA, NA))
qq_results <- rbind(qq_results, c(14, 153851969, NA, NA))
qq_results <- rbind(qq_results, c(15, 157681621, NA, NA))
qq_results <- rbind(qq_results, c(16, 86898991, NA, NA))
qq_results <- rbind(qq_results, c(17, 69701581, NA, NA))
qq_results <- rbind(qq_results, c(18, 61220071, NA, NA))
qq_results <- rbind(qq_results, c(19, 144288218, NA, NA)) #X chrom
qq_results <- rbind(qq_results, c(20, 1637716, NA, NA)) #Y chrom

### Use modified qqman manhattan function to format the data appropriately for a Manhattan plot
d <- manhattan(qq_results, X=19, Y=20)

##Change the chroms from numeric to factors, for colouring purposes.
d$CHR <- as.factor(d$CHR)

#add the end limit to chrY for minor_ticks
minor_ticks <- c(minor_ticks, tail(minor_ticks, 1) + 1637716)

##Add gene names
# geneNames <- complete_results[,c(1,5,22)]
# geneNames$key <- paste(geneNames$snp_id, geneNames$FDR, sep = ":")
# geneNames <- geneNames[,-c(2,3)]
# 
# d$key <- paste(d$SNP, d$P, sep = ":")
# 
# e <- merge(d, geneNames, all = T)

##GGPLOT MOD BY RUSSELL FRASER
#scale_x_continous seems to break the x-axis of the facet_zoom
#There's probably a smart workaround, but for now, the best I can come up with is saving a copy of the graph with
# and without the sacle_x_continous as SVG and then merging the two. 

cols <- c("darkolivegreen4", "dodgerblue3","orange", "blue","orange", "blue","orange", "blue","orange", "blue","orange", "blue","orange", "blue","orange", "blue","orange", "blue","orange", "blue")
cols <- c("darkolivegreen4", "dodgerblue3","darkolivegreen4", "dodgerblue3","darkolivegreen4", "dodgerblue3","darkolivegreen4", "dodgerblue3","darkolivegreen4", "dodgerblue3","darkolivegreen4", "dodgerblue3","darkolivegreen4", "dodgerblue3","darkolivegreen4", "dodgerblue3","darkolivegreen4", "dodgerblue3","darkolivegreen4", "dodgerblue3")


ggplot(d) + 
  geom_point(aes(x=pos, y=logp, color=CHR), size = 1) + 
  theme_bw()+ ylab("-log(p)") + 
  theme(axis.text.x=element_text(angle = 60, hjust=1), legend.title=element_blank(), legend.position = "none") + 
  theme(panel.grid.minor.x = element_line(colour = "light grey"), panel.grid.major.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  # scale_color_viridis(discrete=T) + 
  scale_color_manual(values = cols) +
  geom_hline(yintercept=-log10(0.001), colour="gray33", linetype = 3) + 
  geom_hline(yintercept=-log10(0.05), colour="gray33") + 
  scale_x_continuous(minor_breaks = minor_ticks, breaks = ticks, labels = labs) +
  # scale_x_continuous(labels = scientific_10) +
  # facet_zoom(x=CHR==14 & BP > 20000000 & BP < 110000000) +
  theme(strip.background = element_rect(fill="grey88")) +
  xlab("Chromosome") +
  theme(text = element_text(size = 10)) #+
  # geom_text_repel(data=subset(d, d$logp >= 3), aes(x=pos, y = logp, label = SNP), nudge_x = 0.1)



ggsave("~/Dropbox/figure1.tiff", dpi = 600, units = "mm", width = 190, height = 100)


scientific_10 <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}
