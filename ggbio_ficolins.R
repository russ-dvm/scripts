biocLite("GenomicFeatures")
library(GenomicFeatures)
library(ggbio)
library(biomaRt)
library(VariantAnnotation)

##This is an attempt to create the figure with ggbio that I wanted to create using Gviz (see Gviz_ficolins.R) but couldn't.


#in-house equine genome created according to the instructions at http://bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf
library("BSgenome.Ecaballus.UCSC.equCab2")
genome = "BSgenome.Ecaballus.UCSC.equCab2"

##A variable to contain transcript ids - under the variable transcript name
txname <-  c("ENSECAT00000026503", "ENSECAT00000000831")

##Set the x-axis limits
limits <- GRanges("chr25", IRanges(36817496, 36825451))

##Useful to remember that the txdb object is like a database, and has to be queried in the same way - using SQL like queries. Use select() to get info from the above generated TxDB object. 
# columns(fcns)
# keytypes(fcns)

##Reference sequence
ref_track <- autoplot(genome, which = limits, geom="point")
ref_track <- ref_track + theme(legend.position = "none" )


##Grab transcript data from biomart
tx_data <- makeTxDbFromBiomart(biomart = "ensembl", dataset = "ecaballus_gene_ensembl", transcript_ids = txname)
#Correct the seqlevels
seqlevels(tx_data) <- "chr25"
tx_track <-autoplot(tx_data, which = limits)
tx_track

##Import variants
vcf <- readVcf("~/equine/2014_11_24/8_variants/jun6_ac4/with_fcn1-like/all/jun6.ac4.all.all.vcf", genome = genome)
vcf <- as(vcf[, 1:3], "VRanges")


vcf_plot <- autoplot(vcf, geom = "text", which = limits)


tracks(tx_track, ref_track) + 
  theme_clear() + 
  theme(legend.position = "none")


####Abandoned -- got GViz to work...
