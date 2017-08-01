# source("https://bioconductor.org/biocLite.R")
# biocLite("Gviz")
library("Gviz")
library("BSgenome")
library("biomaRt")
library("VariantAnnotation")

#in-house equine genome created according to the instructions at http://bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf
library("BSgenome.Ecaballus.UCSC.equCab2")

#Gtrack - genome track. Displays the coordinates of the genome.
gtrack <- GenomeAxisTrack()

#idogram track - provides relative position in the chromosome
itrack <- IdeogramTrack(genome = "equCab2", chr = 25)

#sequence track - can provide actual sequence info
strack <- SequenceTrack(Ecaballus, chromosome = 25)

##UCSC genome name and Ensembl name are incompatible for reasons beyond my comprehension. Can access the proper ensembl DB by specifying the dataset to be used with the following command: 
equcab_mart = useMart("ensembl", dataset="ecaballus_gene_ensembl")

#Now make sure to include the "mart" option in the following, wiht mart = equcab_mart
bmt <- BiomartGeneRegionTrack(genome = "EquCab2", stacking = "dense", biomart = equcab_mart, name = "Ficolins", gene = c("ENSECAG00000000436", "ENSECAG00000024620"))

#vcf-track - this relies on the VariantAnnotation package from BioConductor. Access the GRanges for input into Gviz by using the rowRanges function.
vcf <- readVcf("~/equine/2014_11_24/8_variants/jun6_ac4/with_fcn1-like/all/jun6.ac4.all.all.vcf", genome = "equCab2")
vcf_track <- DataTrack(rowRanges(vcf), name = "Variants")


plotTracks(c(itrack, gtrack, bmt, strack, vcf_track), from = 36817496, to = 36825451, shape = "arrow")
s