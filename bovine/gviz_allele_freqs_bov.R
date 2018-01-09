# source("https://bioconductor.org/biocLite.R")
# biocLite("Gviz")
library("Gviz")
library("BSgenome")
library("biomaRt")
library("VariantAnnotation")
library("BSgenome.Btaurus.UCSC.bosTau8")

##Set variables:
genome <- "bosTau8"
all_bams <- "~/bovine/merged_runs/3_realigned/all.bam"

##Gtrack - genome track. Displays the coordinates of the genome.
gtrack <- GenomeAxisTrack()

## Target track - track containing the actual targeted regions
target <- read.table("~/bovine/ref_files/bait.body.txt", sep = "\t")
target$V1 <- gsub("chr", "", target$V1)
targetTrack <- AnnotationTrack(start = target$V2, end = target$V3, chromosome = target$V1, name = "Baits", stacking = "dense")

##Create the transcript track
btau_mart <- useMart("ensembl", dataset="btaurus_gene_ensembl")

##MASP1
# startMasp <- 80546982; endMasp <- 80650824; nameMasp <- "MASP1"; chrMasp <- 1; justMasp <- "left"

## COLEC LOCUS
# start <- 35543220; end <- 35870565; name <- "Collectin Locus"; chr <- 28; just <- "left"

## COLEC11
# start <- 112867044; end <- 112896491; name <- "COLEC11"; chr <- 8; just <- "right"

## FCN1
# start <- 106773026; end <- 106834643; name <- "FCN1"; chr <- 11; just <- "left";

## COLEC12
start <-35629346; end <- 35866133; name <- "COLEC12"; chr <- 24; just <- "right"

##Now make sure to include the "mart" option in the following, with mart = btau_mart
bmtMasp <- BiomartGeneRegionTrack(genome = "BosTau8", biomart = btau_mart, name = nameMasp, collapseTranscripts = T, chromosome = chrMasp, start = startMasp, end = endMasp, shape = "arrow", transcriptAnnotation = "symbol", just.group = justMasp)

bmt <- BiomartGeneRegionTrack(genome = "BosTau8", biomart = btau_mart, name = name, collapseTranscripts = T, chromosome = chr, start = start, end = end, shape = "arrow", transcriptAnnotation = "symbol", just.group = just)
##Need to fix the bmt track, as the coordinates for certain genes are not the same as what we used (eg the collectin locus -sfpa/mbl)
## Convert to data frame and get rid of the fluff
# IF FCN1!!!: 
# symbol(bmt)[c(1:21)] <- "FCN1"

geneDf <- as.data.frame(ranges(bmt))
#Remove fake-out SPA
geneDf <- geneDf[-grep("ENSBTAT00000064451", geneDf$transcript),]
geneDf <- geneDf[-grep("ENSBTAT00000063315", geneDf$transcript),]
#Remove fake-out CGN
geneDf <- geneDf[-grep("ENSBTAT00000025788", geneDf$transcript),]
#Label genes
geneDf[grep("ENSBTAT00000001165", geneDf$transcript),]$symbol <- "MBL1"
geneDf[grep("ENSBTAG00000006536", geneDf$gene),]$symbol <- "CGN1"
geneDf[grep("ENSBTAG00000046421", geneDf$gene),]$symbol <- "SFTPD"
#Assign MBL1 a new gene id
geneDf[grep("ENSBTAT00000001165", geneDf$transcript),]$gene <- "MBL1"
#Fix funky CL43 coords
geneDf[grep("ENSBTAG00000047317", geneDf$gene),]$symbol <- "CL43"
geneDf[grep("ENSBTAE00000030432", geneDf$exon),]$start <- 35722932
geneDf[grep("ENSBTAE00000030432", geneDf$exon),]$end <- 35723104


geneTrack <- GeneRegionTrack(geneDf, genome = genome, collapseTranscripts = T, shape = "arrow", transcriptAnnotation = "symbol", name = "", chromosome = chr)



#### ADD AN ALLELE FREQUENCY TRACK ####
freqs <- read.table("~/Dropbox/temp/pub_table.txt", h=T, sep="\t")
##Slim down the table
freqsSlim <- freqs[,c(1,2,13)]
freqFunc <- function(x, start, end, chr) {
  freqsSlimGene <- subset(freqsSlim, Position >= start & Position <= end & Chrom == paste("chr", chr, sep = ""))
  #TBC after iRanges
  
  # Create the IRanges
  irangeGene <- c(start:end)
  width <- rep(1, length(irangeGene))
  freqIrange <- IRanges(start = irangeGene, end = irangeGene, width = width)
  seqnames <- chr
  
  ## Finish the freq table
  tmpDF <- data.frame("Position" = irangeGene, "Chrom" = paste("chr", chr, sep = ""), "logBH" = NA)
  tmpMerged <- merge(tmpDF, freqsSlimGene, by = "Position", all.x = T)
  tmpMerged <- tmpMerged[,c(1,2,5)]
  colnames(tmpMerged) <- c("Position", "Chrom", "logBH")
  
  ####Build the GRanges object####
  freqGrange <- GRanges(seqnames = seqnames, ranges = freqIrange, strand = NULL)
  mcols(freqGrange)$logBH <- tmpMerged$logBH
  
  
  ##Make the datatrack
  class2 <- -log10(5e-6)
  
  DataTrack(freqGrange, type = "p", name = "-log(p-adj)", chromosome = chr, showSampleNames = T, baseline = class2, col.baseline = "red", size = 1, grid = F, ylim = c(0,8), fontsize = 11)
}

maspVar <- freqFunc(freqsSlim, startMasp, endMasp, chrMasp)
colecVar <- freqFunc(freqsSlim, start, end, chr)

##Depth track - to allow identification of low coverage regions where variants may not have been called
dtrack <- DataTrack(range = all_bams, genome = genome, type = "histogram", col.histogram = "#377EB8", fill = "#377EB8", name = "Depth", baseline = 600, col.baseline = "red", size = 1, chromosome = chr)

dtrackMasp <- DataTrack(range = all_bams, genome = genome, type = "histogram", col.histogram = "#377EB8", fill = "#377EB8", name = "Depth", baseline = 600, col.baseline = "red", size = 1, chromosome = chrMasp)


#The boxes to the left are the ".title" boxes - thus to manipulate the colours, use things like background.title or fontcolor.title. To see the names of all the parameters that can be manipulated -- names(displayPars(**TRACK**)). Can manipulate globally (within the plotTracks), or for individual tracks (within the tracks themselves). background.title = "white", fontcolor.title = "black"

## Manipualte the plotting device to create a composite.
tiff("~/Dropbox/temp/test.tiff", width = 6.85, height = 6, units = "in", res = 1000)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,1)))

pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
plotTracks(c(gtrack, bmtMasp, maspVar, dtrackMasp, targetTrack), from = startMasp, to = endMasp, sizes = c(0.5,0.5,1,1,0.5), background.title = "white", fontcolor.title = "white", col.axis = "black", fontsize = 10, cex.axis = 0.7, cex.title = 0.7, add = T)

popViewport(1)
pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))

plotTracks(c(gtrack, geneTrack, colecVar, dtrackMasp, targetTrack), from = start, to = end, sizes = c(0.5,0.5,1,1,0.5), background.title = "white", fontcolor.title = "white", col.axis = "black", fontsize = 10, cex.axis = 0.7, cex.title = 0.7, add = T)

dev.off()

## For supp figs:
plotTracks(c(gtrack, bmt, colecVar, dtrack, targetTrack), from = start, to = end, sizes = c(0.5,0.5,1,1,0.5), background.title = "white", fontcolor.title = "black", col.axis = "black", fontsize = 10, cex.axis = 0.7, cex.title = 0.7, add = T)
