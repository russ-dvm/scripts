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

# # ##MASP1
start <- 80546982
end <- 80650824
name <- "MASP1"
chr <- 1
just <- "left"

## COLEC LOCUS
# start <- 35543220
# end <- 35870565
# name <- "Collectin Locus"
# chr <- 28
# just <- "left"

# ## COLEC11
# start <- 112867044
# end <- 112896491
# name <- "COLEC11"
# chr <- 8
# just <- "right"

# ## FCN1
# start <- 106773026
# end <- 106834643
# name <- "FCN1"
# chr <- 11
# just <- "left"

# ## COLEC12
# start <-35629346		
# end <- 35866133
# name <- "COLEC12"
# chr <- 24
# just <- "right"

##Gtrack - genome track. Displays the coordinates of the genome.
gtrack <- GenomeAxisTrack()

## Target track - track containing the actual targeted regions
target <- read.table("~/bovine/ref_files/bait.body.txt", sep = "\t")
target$V1 <- gsub("chr", "", target$V1)
targetTrack <- AnnotationTrack(start = target$V2, end = target$V3, chromosome = target$V1, name = "d", stacking = "dense")

##Depth track - to allow identification of low coverage regions where variants may not have been called
dtrack <- DataTrack(range = all_bams, chromosome = chr, genome = genome, type = "histogram", col.histogram = "#377EB8", fill = "#377EB8", name = "Depth", baseline = 600, col.baseline = "red", size = 1)

##Create the transcript track
##UCSC genome name and Ensembl name are incompatible for reasons beyond my comprehension. Can access the proper ensembl DB by specifying the dataset to be used with the following command: 
btau_mart <- useMart("ensembl", dataset="btaurus_gene_ensembl")

##Now make sure to include the "mart" option in the following, with mart = btau_mart;  stacking = "squish"
bmt <- BiomartGeneRegionTrack(genome = "BosTau8", biomart = btau_mart, name = "", collapseTranscripts = T, chromosome = chr, start = start, end = end, shape = "arrow", transcriptAnnotation = "symbol", just.group = just)

##Need to fix the bmt track, as the coordinates for certain genes are not the same as what we used (eg the collectin locus -sfpa/mbl)
## Convert to data frame and get rid of the fluff

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

var_track <- DataTrack(freqGrange, type = "p", name = "-log(p-adj)", chromosome = chr, showSampleNames = T, baseline = class2, col.baseline = "red", size = 1, grid = F, ylim = c(0,8), fontsize = 11)

#The boxes to the left are the ".title" boxes - thus to manipulate the colours, use things like background.title or fontcolor.title. To see the names of all the parameters that can be manipulated -- names(displayPars(**TRACK**)). Can manipulate globally (within the plotTracks), or for individual tracks (within the tracks themselves). background.title = "white", fontcolor.title = "black"

# plotTracks(c(gtrack, geneTrack, var_track, dtrack), from = start, to = end, sizes = c(0.5,0.5,1,0.5), background.title = "white", fontcolor.title = "black", col.axis = "black", fontsize = 10, cex.axis = 0.7, cex.title = 0.7)

plotTracks(c(gtrack, bmt, var_track, dtrack, targetTrack), from = start, to = end, sizes = c(0.5,0.5,1,1,0.5), background.title = "white", fontcolor.title = "black", col.axis = "black", fontsize = 10, cex.axis = 0.7, cex.title = 0.7)

## IF FCN1!!!: symbol(bmt)[c(1:21)] <- "FCN1"

