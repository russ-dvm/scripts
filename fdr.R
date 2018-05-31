## Function that performes both the base R p.adjust, as well as manually calculating the BH critical value.

fdr <-
function(x, col, fdr){
  #determine which column contains the P values
  pcol <- which(colnames(x) == col)

  #total number of tests
  m <- nrow(x)

  #order the columns by increasing p value
  x.ordered <- x[order(x[,pcol]),]

  #assign ranks, smallest pvalue = 1
  x.ordered$bh.rank<- 1:m

  #determine the critical value (iQ/m) (ie rank*fdr/total)
  x.ordered$bh.crit.val <- x.ordered$bh.rank*fdr/m

  #a p value is significant if it is less than the critical value
  x.ordered$bh.sig <- x.ordered[,pcol] < x.ordered$bh.crit.val

  #In fact, all p values ranked LOWER (ie a smaller number) than the highest p value that is < the critical value are considered significant.
  # ie if the 10th pvalue < critical value, then pvalues1-10 are ALL sig.
  k <-  max(x.ordered[x.ordered$bh.sig == T,]$bh.rank)
  x.ordered$bh.sig.final <- ifelse(x.ordered$bh.rank <= k, T, F)
  
  #perform the base R p.adjust. Here FDR is used to determine the significance cutoff.
  x.ordered$p.adj <- p.adjust(x.ordered[,pcol], method = "BH")
  x.ordered$p.adj.sig <- x.ordered$p.adj < fdr

  x.ordered$method.comp <- x.ordered$bh.sig.final == x.ordered$p.adj.sig

  x.ordered
}
