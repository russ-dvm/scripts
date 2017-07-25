#RS Fraser, 2017-07-24
#Check the correlation of probes for the same gene from Hein's microarray.

library(tidyverse)
library(psych)


innate <- read.table("~/Desktop/innate_immune_gene_expression_NOT_GAPDH_norm_for_correlation.txt", h=T)


#Make a susbet containing only probes that have a duplicate
dup_probes <- innate[,grep(".p", colnames(innate))]

#Further subset so that only a single probe for each gene is included
dup_probe_1 <- dup_probes[,grep("*.p1", colnames(dup_probes))]

#Create an empty list for the results
results <- list()

#Loop over the data frame. 
for (i in 1:ncol(dup_probe_1)){
  
  #Determine gene name. Note that occasionally there will be overlap between two genes (e.g. C5 and C5AR1). Although technically incorrect, it's easy enough to look past this.
  gene_name <- gsub(".p1", "", colnames(dup_probe_1)[i])
  
  #Make sure the loop goes through the data frame with ALL the probes
  dups <- grep(gene_name, colnames(dup_probes))
  ct <- corr.test(dup_probes[,dups])
  a <- ct$r
  results[[i]] <- a
  
}

results


##Check the heatmap values
bd129 <- innate[,grep("BD129", colnames(innate))]
bd129$avg <- rowMeans(bd129)

c1s <- innate[,grep("C1S", colnames(innate))]
c1s$avg <- rowMeans(c1s)

clec7a <- innate[, grep("CLEC7A", colnames(innate))]
clec7a$avg <- rowMeans(clec7a)

tlr3 <- innate[,grep("TLR3", colnames(innate))]
tlr3


#Values do not match Hein's.... could they be GAPDH normalized?
gapdh <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/gapdh.txt", h=T)
gapdh <- melt(gapdh)
gapdh$id <- innate$id
bd129$id <- innate$id
bd129_gap <- merge(bd129, gapdh, by = "id")
bd129_gap$p1_gap <- bd129_gap$BD129.p1 - bd129_gap$value
bd129_gap$p2_gap <- bd129_gap$BD129.p2 - bd129_gap$value
bd129_gap$gap_avg <- rowMeans(bd129_gap[,8:9])
#Nope, not due to log transformation.

##BL figured out that the averages were performed as follows: the raw data (not available here) was averaged and THEN log2 transformed. The values used here are all log2 transformed; thus 
bd129_gap$log <- log(bd129_gap$avg)



##BL requested check of ref genes and acute phase protein genes, as well as a comparison of APP genes vs top 20 variable innate immune genes
heatmap <- read.table("~/Dropbox/Research/Lab Book/NGS/oink/microarray/figs_for_paper/heatmap_data.txt", h=T, sep = "\t")

#APP & Ref genes
app <- heatmap[,c("CP.AVG", "CRP", "HP.AVG", "ITIH4.AVG", "ORM1", "SAA.AVG")]
corr.test(app)

ref <- heatmap[, c("ACTB.AVG", "B2M.AVG", "GAPDH", "HPRT1", "SDHA.AVG")]
corr.test(ref)

#Comparison of ITIH4 and SAA to top 20 innate immune genes
top20 <- c("DDX3Y", "MBL2", "SCGB1A1", "PR39", "PMAP.23", "PMAP.37", "NPG4.AVG", "SFTPD", "PGLYRP1.A", "HAMP.AVG", "PGLYRP2.B", "MYD88", "CLEC1B", "LBP", "PMAP.36", "PBD.2", "ITLN2.AVG", "KLRF1", "DDX58.AVG","BD.104.like.2")

comp_top_20 <- function(app, top20) {
  app_list <- list()
  for (i in 1:length(top20)){
    a <- heatmap[, c(app, top20[i])]
    b <- corr.test(a)
    c <- b$r
    app_list[[i]] <- c
  }
  return(app_list)
}

itih_list <- comp_top_20("ITIH4.AVG", top20)
saa_list <- comp_top_20("SAA.AVG", top20)
cp_list <- comp_top_20("CP.AVG", top20)
hp_list <- comp_top_20("HP.AVG", top20)
crp_list <- comp_top_20("CRP", top20)
orm_list <- comp_top_20("ORM1", top20)




#some graphs for BL
ggplot(heatmap, aes(x = ITIH4.AVG, y = PGLYRP2.B)) + geom_point() + geom_smooth(method = lm, se = T)
ggplot(heatmap, aes(x = ITIH4.AVG, y = MYD88)) + geom_point() + geom_smooth(method = lm, se = T)
ggplot(heatmap, aes(x = ITIH4.AVG, y = LBP)) + geom_point() + geom_smooth(method = lm, se =T)
ggplot(heatmap, aes(x = ITIH4.AVG, y = DDX58.AVG)) + geom_point() + geom_smooth(method = lm, se = T)
ggplot(heatmap, aes(x = SAA.AVG, y = LBP)) + geom_point() + geom_smooth(method = lm, se = T)
