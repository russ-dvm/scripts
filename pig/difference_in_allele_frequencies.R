library(tidyverse)

viro <- read.table("~/Dropbox/temp/viro.tsv", h = T, sep = "\t")
bact <- read.table("~/Dropbox/temp/bact.tsv", h = T, sep = "\t")
cond <- read.table("~/Dropbox/temp/cond.tsv", h = T, sep = "\t")

ratios <- function(x) {
  colnames(x) <- gsub("\\.*n*\\d*", "", colnames(x))
  for (i in 2:ncol(x)) {
    a <- colnames(x)[i]
    b <- paste("rs.", a, sep = "")
    x[,b] <- x[,i]/x$Healthy
  }
  y <- x[,grep("rs", colnames(x))]
  y <- y[,-2]
  z <- gather(y, rsID)
  colnames(z) <- c("rsID", "comp", "value")
  z[,2] <- gsub("rs.", "", z[,2])
  return(z)
}

viroR <- ratios(viro)
bactR <- ratios(bact)
condR <- ratios(cond)


ggplot(viroR, aes(x = rsID, y = log2(value))) + theme_classic() + geom_boxplot() + geom_point(aes(colour = comp)) + theme(axis.text.x = element_text(hjust = 1, angle = 60)) + theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_rect(linetype = "solid", colour = "black")) + geom_hline(yintercept = 0, linetype = "dashed")
ggplot(bactR, aes(x = rsID, y = log2(value))) + theme_classic() + geom_boxplot() + geom_point(aes(colour = comp)) + theme(axis.text.x = element_text(hjust = 1, angle = 60)) + geom_hline(yintercept = 0, linetype = "dashed")
ggplot(condR, aes(x = rsID, y = log2(value))) + theme_classic() + geom_boxplot() + geom_point(aes(colour = comp)) + theme(axis.text.x = element_text(hjust = 1, angle = 60)) + geom_hline(yintercept = 0, linetype = "dashed")

