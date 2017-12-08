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

# ##MASP1
# start <- 80546982
# end <- 80650824
# name <- "MASP1"
# chr <- 1
# geneId <- "ENSBTAG00000012467"

## COLEC LOCUS
start <- 35543220
end <- 35870565
name <- "Collectin Locus"
chr <- 28
transID <- c("ENSBTAT00000018649", "ENSBTAT00000028716", "ENSBTAT00000003773", "ENSBTAT00000008579", "ENSBTAT00000001165", "ENSBTAT00000031298")

##Gtrack - genome track. Displays the coordinates of the genome.
gtrack <- GenomeAxisTrack()

##ideogram track - provides relative position in the chromosome
# itrack <- IdeogramTrack(genome = genome, chromosome = chr)

##Depth track - to allow identification of low coverage regions where variants may not have been called
dtrack <- DataTrack(range = all_bams, chromosome = chr, genome = genome, type = "histogram", col.histogram = "#377EB8", fill = "#377EB8", name = "Depth of Coverage (all samples)", baseline = 600, col.baseline = "red", size = 1)

##Create the transcript track
##UCSC genome name and Ensembl name are incompatible for reasons beyond my comprehension. Can access the proper ensembl DB by specifying the dataset to be used with the following command: 
btau_mart = useMart("ensembl", dataset="btaurus_gene_ensembl")

##Now make sure to include the "mart" option in the following, with mart = equcab_mart
bmt <- BiomartGeneRegionTrack(genome = "BosTau8", stacking = "squish", biomart = btau_mart, name = name, collapseTranscript = T, chromosome = chr, start = start, end = end)

plotTracks(bmt)

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
class1 <- -log10(2e-7)
class2 <- -log10(5e-6)

var_track <- DataTrack(freqGrange, name = "Allele Frequency", type = "p", chromosome = chr, showSampleNames = T, baseline = class2, col.baseline = "red", size = 1, grid = F, ylim = c(0,8))

##Highlight track - highlight exons? regions of increased variation? Note that the tracks included in the highlight track will need to be removed from the final plotTrack call.
#htrack <- HighlightTrack(trackList = list(bmt, dtrack, vcf_track), start = 36817496, width = 1000, chromosome = chr)


#The boxes to the left are the ".title" boxes - thus to manipulate the colours, use things like background.title or fontcolor.title. To see the names of all the parameters that can be manipulated -- names(displayPars(**TRACK**)). Can manipulate globally (within the plotTracks), or for individual tracks (within the tracks themselves).
plotTracks(c(gtrack, bmt, var_track, dtrack), from = start, to = end, min.height = 5)
plotTracks(c(gtrack, bmt, var_track, itrack), from = start, to = end, min.height = 5)


