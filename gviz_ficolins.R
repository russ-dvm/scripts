# source("https://bioconductor.org/biocLite.R")
# biocLite("Gviz")
library("Gviz")
library("BSgenome")
library("biomaRt")
library("VariantAnnotation")

#in-house equine genome created according to the instructions at http://bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf
library("BSgenome.Ecaballus.UCSC.equCab2")

##Set variables:
genome <- "equCab2"
all_bams <- "~/equine/2014_11_24/6_realigned/all.realigned.bam"

##FCN-1
start <- 36817400
end <- 36825451
name <- "FCN1"
chr <- 25

##FCN-1-like
# start <- 36801475
# end <- 36804532
# name <- "FCN-1-like"
# chr <- 25

##Gtrack - genome track. Displays the coordinates of the genome.
gtrack <- GenomeAxisTrack()

##ideogram track - provides relative position in the chromosome
itrack <- IdeogramTrack(genome = genome, chromosome = chr)

##sequence track - can provide actual sequence info
strack <- SequenceTrack(Ecaballus, chromosome = chr)

##Depth track - to allow identification of low coverage regions where variants may not have been called
dtrack <- DataTrack(range = all_bams, chromosome = chr, genome = genome, type = "histogram", col.histogram = "#377EB8", fill = "#377EB8", name = "Depth of Coverage (all samples)", baseline = 890, col.baseline = "red", size = 1)

##Create the transcript track
##UCSC genome name and Ensembl name are incompatible for reasons beyond my comprehension. Can access the proper ensembl DB by specifying the dataset to be used with the following command: 
equcab_mart = useMart("ensembl", dataset="ecaballus_gene_ensembl")

##Now make sure to include the "mart" option in the following, with mart = equcab_mart
bmt <- BiomartGeneRegionTrack(genome = "EquCab2", stacking = "squish", biomart = equcab_mart, name = name, gene = c("ENSECAG00000000436", "ENSECAG00000024620"))

##The "histogram" track really appears to be more of a bargraph to me. I cannot get it to plot frequencies/counts, just absolute values... so not really sure. I did post in the Bioconductors forum about this with no response. https://support.bioconductor.org/p/98774/ 
##Decided to try an alternate strategy: use the annotated_variant dataset generated for "rate_of_snps.R" and see if I can summarize and then convert it into a GRanges object, so that the bargraph/"histogram" will actually represent what we want.

####SIMPLE graph - plots the total number of variants without any annotation
####COMPLICATED graph - colours the histogram according to the different variant types. Note that because SnpEff annotates based on the information for multiple genes, the annotation with respect to OUR gene of interest is sometimes incorrect. E.g. for the Ficolins, you may have a SNP labeled as "downstream" when it is in fact an intron snp - this is because the snp occurs downstream of a gene that is within 5kb of the FCN...Anyways, can use the "feature" info of the annotated_depth DF to correct this, I think. Use the SnpEff info for better annotation on non-synonymous codingifelse variant...
##Simple
annotated_depth <- read.table("~/equine/2014_11_24/depth_of_regions/results_utr_adjusted.txt", h=T, sep="\t")
annotated_depth <- read.table("~/Desktop/annotated_depth_and_variants.txt", h=T, sep="\t")
##Keep only variants
annotated_depth_variants <- annotated_depth[annotated_depth$is_variant == T,]
##Slim down the table
variants <- annotated_depth_variants[,c(1:2,5,12)]
##Only the variants in the selected gene
gene_variants <- subset(variants, pos >= start & pos <= end)
gene_missense_variants <- gene_variants[grep("missense", gene_variants$variant_type),]
gene_stop_variants <- gene_variants[grep("stop", gene_variants$variant_type),]

##Set bin width for the histogram function
bin_width <- 50

## Pass info from the hist function to a variable
histinfo <- hist(gene_variants$pos, breaks = seq(from = start, to = end, by = bin_width))
histinfo_missense <- hist(gene_missense_variants$pos, seq(from = start, to = end, by = bin_width))
histinfo_stop <- hist(gene_stop_variants$pos, seq(from = start, to = end, by = bin_width))
breaks <- histinfo$breaks
var_counts <- histinfo$counts
var_counts_mis <- histinfo_missense$counts
var_counts_stop <- histinfo_stop$counts
#Breaks are boundaries, so there is usually the #breaks - 1 number of counts. Need to add an extra empty count to make Gviz happy downstream.
var_counts <- c(var_counts, 0)
var_counts_mis <- c(var_counts_mis, 0)
var_counts_stop <- c(var_counts_stop, 0)
var_counts_not_mis <- var_counts - var_counts_mis - var_counts_stop

##Complicated
# ann_var <- read.table("~/equine/2014_11_24/depth_of_regions/temp2.txt", h=T, sep="\t")
# ann_var <- ann_var[,c(1:2,12)]
# ann_var_gene <- subset(ann_var, pos >= start & pos <= end)

###COMPLICATED - not toally working. Issues with the mutliple annotation types.
# bin_width <- 50
# 
# hist_subfunction <- function(x, data, start, end, bin_width) {
#   a <- subset(data, data$variant_type == x)
#   b <- hist(a$pos, breaks = seq(from = start, to = end, by = bin_width))
#   return(b)
# }
# 
# var_types <- as.vector(unique(ann_var_gene$variant_type))
# list_of_hist_info <- sapply(var_types, hist_subfunction, data = ann_var_gene, start = start, end = end, bin_width = bin_width, simplify = F, USE.NAMES = T)
# list_of_hist_info$synonymous_variant$counts <- c(list_of_hist_info$synonymous_variant$counts, 0)
# list_of_hist_info$missense_variant$counts <- c(list_of_hist_info$missense_variant$counts, 0)
# list_of_hist_info$`splice_region_variant&intron_variant`$counts <- c(list_of_hist_info$`splice_region_variant&intron_variant`$counts, 0)
# list_of_hist_info$downstream_gene_variant$counts <- c(list_of_hist_info$downstream_gene_variant$counts, 0)
# list_of_hist_info$stop_gained$counts <- c(list_of_hist_info$stop_gained$counts, 0)
# list_of_hist_info$`missense_variant&splice_region_variant`$counts <- c(list_of_hist_info$`missense_variant&splice_region_variant`$counts, 0)

##############################
###Build the GRanges object###
##############################

##Start with the IRange
start_irange <- breaks
end_irange <- start_irange + (bin_width - 1)
width <- rep(bin_width, length(start_irange))

variant_irange <- IRanges(start = start_irange, end = end_irange, width = width)
seqnames <- unique(variants$chrom)
seqname_fcn <- "chr25"

var_grange <- GRanges(seqnames = seqname_fcn, ranges = variant_irange, strand = NULL)

##SIMPLE GRAPH - this adds only the total variant counts to the figure.
mcols(var_grange)$var_counts <- var_counts

##IN BETWEEN - separates missense and non-missense.
mcols(var_grange)$synonymous_or_intron <- var_counts_not_mis
mcols(var_grange)$missense <- var_counts_mis
mcols(var_grange)$stop <- var_counts_stop

##ADDED ANNOTATIONS - order is important - keep order in the DataTrack below.
# mcols(var_grange)$synonymous <- list_of_hist_info$synonymous_variant$counts
# mcols(var_grange)$downstream <- list_of_hist_info$downstream_gene_variant$counts
# mcols(var_grange)$missense <- list_of_hist_info$missense_variant$counts
# mcols(var_grange)$missense_splice <- list_of_hist_info$`missense_variant&splice_region_variant`$counts
# mcols(var_grange)$stop <- list_of_hist_info$stop_gained$counts
# mcols(var_grange)$intron <- list_of_hist_info$`splice_region_variant&intron_variant`$counts
  

##Make the datatrack
var_track <- DataTrack(var_grange, name = "Variants", type = "histogram", chromosome = chr, showSampleNames = T, groups = c("other", "missense", "stop"))

##Highlight track - highlight exons? regions of increased variation? Note that the tracks included in the highlight track will need to be removed from the final plotTrack call.
#htrack <- HighlightTrack(trackList = list(bmt, dtrack, vcf_track), start = 36817496, width = 1000, chromosome = chr)

#The boxes to the left are the ".title" boxes - thus to manipulate the colours, use things like background.title or fontcolor.title. To see the names of all the parameters that can be manipulated -- names(displayPars(**TRACK**)). Can manipulate globally (within the plotTracks), or for individual tracks (within the tracks themselves).
plotTracks(c(itrack, gtrack, bmt, var_track, dtrack, strack), from = start, to = end, min.height = 5)
plotTracks(c(itrack, gtrack, bmt, var_track, strack), from = start, to = end, min.height = 5)


