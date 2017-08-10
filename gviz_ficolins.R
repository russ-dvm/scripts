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
# start <- 36817496
# end <- 36825451
# name <- "FCN1"
# chr <- 25

##FCN-1-like
start <- 36801475
end <- 36804532
name <- "FCN-1-like"
chr <- 25

##Gtrack - genome track. Displays the coordinates of the genome.
gtrack <- GenomeAxisTrack()

##ideogram track - provides relative position in the chromosome
itrack <- IdeogramTrack(genome = genome, chromosome = chr)

##sequence track - can provide actual sequence info
strack <- SequenceTrack(Ecaballus, chromosome = chr)

##Depth track - for fun
dtrack <- DataTrack(range = all_bams, chromosome = chr, genome = genome, type = "histogram", col.histogram = "#377EB8", fill = "#377EB8", name = "Depth of Coverage (all samples)", baseline = 890, col.baseline = "red", size = 1)

#altrack <- AlignmentsTrack(range = all_bams, chromosome = chr, genome = genome, name = "alignment")

##UCSC genome name and Ensembl name are incompatible for reasons beyond my comprehension. Can access the proper ensembl DB by specifying the dataset to be used with the following command: 
equcab_mart = useMart("ensembl", dataset="ecaballus_gene_ensembl")

##Now make sure to include the "mart" option in the following, with mart = equcab_mart

bmt <- BiomartGeneRegionTrack(genome = "EquCab2", stacking = "dense", biomart = equcab_mart, name = name, gene = c("ENSECAG00000000436", "ENSECAG00000024620"), shape = "arrow")

##vcf-track - this relies on the VariantAnnotation package from BioConductor. Access the GRanges for input into Gviz by using the rowRanges function.
#vcf <- readVcf("~/equine/2014_11_24/8_variants/jun6_ac4/with_fcn1-like/all/jun6.ac4.all.all.vcf", genome = genome)

##Window = the NUMBER OF BINS, NOT the SIZE of bins... so a lower window value will result in bigger bins.
#vcf_track <- DataTrack(rowRanges(vcf), name = "Variants", type = "histogram", chromosome = chr, window = 200)
##The "histogram" track really appears to be more of a bargraph to me. I cannot get it to plot frequencies/counts, just absolute values... so not really sure. I did post in the Bioconductors forum about this with no response. https://support.bioconductor.org/p/98774/ 

##Decided to try an alternate strategy: use the annotated_variant dataset generated for "rate_of_snps.R" and see if I can summarize and then convert it into a GRanges object, so that the bargraph/"histogram" will actually represent what we want.
annotated_depth <- read.table("~/equine/2014_11_24/depth_of_regions/annotated_depth_and_variants.txt", h=T, sep="\t")
##Keep only variants
annotated_depth_variants <- annotated_depth[annotated_depth$is_variant == T,]
##Slim down the table
variants <- annotated_depth_variants[,1:2]
#Only the ficolin variants
fcn_variants <- subset(variants, pos >= start & pos <= end)
#Get the binning info from the hist function
histinfo <- hist(fcn_variants$pos, breaks = 200)
bin_width <- histinfo$breaks[2] - histinfo$breaks[1]
breaks <- histinfo$breaks
var_counts <- histinfo$counts

min_break <- min(breaks)
max_break <- max(breaks)
break_num <- length(breaks)

#The breaks are only constructed around areas in which there are values - so if start/stop (defined above) have no variants then bins will not be created -this leaves gaps in the final track

##Add the missing bins and data (counts will be 0)

if (end - max_break >= 0) {
  end_difference <- end - max_break
  no_missing_end_bins <- ceiling(end_difference/bin_width)
  missing_end_breaks <- seq((max_break+50), by = bin_width, length.out = no_missing_end_bins)
  breaks <- c(breaks, missing_end_breaks)
  
  ##adjust the counts, too...
  var_counts <- c(var_counts, rep(0, no_missing_end_bins + 1))
}

if (min_break - start >= 0) {
  start_difference <- min_break - start
  no_missing_start_bins <- ceiling(start_difference/bin_width)
  missing_start_breaks <- seq((min_break-50), by = -bin_width, length.out = no_missing_start_bins)
  breaks <- c(rev(missing_start_breaks), breaks)
  
  ##adjust the counts, too...
  var_counts <- c(rep(0, no_missing_start_bins), var_counts)
}



###Build the GRanges object
##Start with the IRange
start_irange <- breaks
end_irange <- start_irange + (bin_width - 1)
width <- rep(bin_width, length(start_irange))

variant_irange <- IRanges(start = start_irange, end = end_irange, width = width)
seqnames <- unique(variants$chrom)
seqname_fcn <- "chr25"

var_grange <- GRanges(seqnames = seqname_fcn, ranges = variant_irange, strand = NULL)
mcols(var_grange)$var_counts <- var_counts

##Make the datatrack
var_track <- DataTrack(var_grange, name = "Variants", type = "histogram", chromosome = chr)

##Highlight track - highlight exons? regions of increased variation? Note that the tracks included in the highlight track will need to be removed from the final plotTrack call.
#htrack <- HighlightTrack(trackList = list(bmt, dtrack, vcf_track), start = 36817496, width = 1000, chromosome = chr)

#The boxes to the left are the ".title" boxes - thus to manipulate the colours, use things like background.title or fontcolor.title. To see the names of all the parameters that can be manipulated -- names(displayPars(**TRACK**)). Can manipulate globally (within the plotTracks), or for individual tracks (within the tracks themselves).
plotTracks(c(itrack, gtrack, bmt, var_track, dtrack, strack), from = start, to = end)

