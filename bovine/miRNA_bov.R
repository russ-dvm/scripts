library(tidyverse)

mir <- read.table("~/bovine/merged_runs/miR_predictions/3kb_from_stop/intersection/intersection.same-miR.txt", h = F, sep = "\t")
mir <- mir[,c(1:5)]
colnames(mir) <- c("chr", "start", "end", "fasta", "miR")
mir$chr <- as.factor(mir$chr)

mirVar <- read.table("~/bovine/merged_runs/miR_predictions/3kb_from_stop/intersection/snps.intersection.vcf", h = F, sep = "\t")
mirVar <- mirVar[,c(1:5, 35:38)]
colnames(mirVar) <- c("chrom", "pos", "rsid", "ref", "alt", "start", "stop", "fasta", "mir")

geneify <- function(x, fasta){
  x$gene <- NA
  try(x[grep("chromosome:UMD3.1:11:106831644:106834644:1", x$fasta),]$gene <- "FCN1")
  try(x[grep("chromosome:UMD3.1:8:112867044:112870044:-1", x$fasta),]$gene <- "COLEC11")
  try(x[grep("chromosome:UMD3.1:26:6348278:6351278:1", x$fasta),]$gene <- "MBL2")
  try(x[grep("chromosome:UMD3.1:28:35722932:35725932:1", x$fasta),]$gene <- "CL43")
  try(x[grep("chromosome:UMD3.1:28:35837918:35840918:-1", x$fasta),]$gene <- "MBL1")
  try(x[grep("chromosome:UMD3.1:24:35629346:35632346:-1", x$fasta),]$gene <- "COLEC12")
  try(x[grep("chromosome:UMD3.1:28:35700096:35703096:1", x$fasta),]$gene <- "CL46")
  try(x[grep("chromosome:UMD3.1:28:35602674:35605674:1", x$fasta),]$gene <- "CGN1")
  try(x[grep("chromosome:UMD3.1:28:35824384:35827384:1", x$fasta),]$gene <- "SFTPD")
  try(x[grep("chromosome:UMD3.1:14:47260661:47263661:-1", x$fasta),]$gene <- "COLEC10")
  try(x[grep("chromosome:UMD3.1:28:35848716:35851716:-1", x$fasta),]$gene <- "SFTPA1")
  try(x[grep("chromosome:UMD3.1:16:43477400:43480400:1", x$fasta),]$gene <- "MASP2")
  try(x[grep("chromosome:UMD3.1:1:80647824:80650824:1", x$fasta),]$gene <- "MASP1")
  return(x)
}

mir <- geneify(x = mir, fasta = fasta)
mirVar <- geneify(x = mirVar, fasta = fasta)

## Run rate_of_snps_bov.R to obtain the colecTrimmed DF.
totSNP <- subset(colecTrimmed, Region == "Downstream 3 kb" & Gene != "Total")
mircount <- data.frame(table(mirVar$gene))
totSNP <- merge(totSNP, mircount, by.x = "Gene", by.y = "Var1", all = T)

mirVar$gene <- factor(mirVar$gene, levels = mircount[order(mircount$Freq, decreasing = T), ]$Var1)

ggplot(mir, aes(x = gene)) + geom_histogram(stat="count")
ggplot(mirVar, aes(x = gene)) + geom_histogram(stat = "count") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) + 
  ylab("Number of MREs containing SNVs") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  xlab("")

ggplot(totSNP, aes(x = Gene, y = Freq/No.variants)) + geom_bar(stat = "identity")
