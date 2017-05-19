

strack <- SequenceTrack(Ecaballus, chromosome = 3)
plotTracks(strack, from = 88947828, to = 88947868, cex = 1, add53 = T, noLetters = T)


from <- 88947828
to <- 89020284

knownGenes <- UcscTrack(genome = "equCab2", chromosome = "chr1", track = "refGene", from = from, to = to, trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", gene = "name", symbol = "name", transcript = "name", strand = "strand", fill = "red", name = "Horse Genes")
otherGenes <- UcscTrack(genome = "equCab2", chromosome = "chr1", track = "Other RefSeq", from = from, to = to, trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", gene = "name", symbol = "name", transcript = "name", strand = "strand", fill = "red", name = "Other Genes")
otherGenes <- UcscTrack(genome = "equCab2", chromosome = "chr1", track = "Human Proteins", from = from, to = to, trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", gene = "name", symbol = "name", transcript = "name", strand = "strand", fill = "red", name = "UCSC Genes")



alTrack <- AlignmentsTrack(system.file(package = "Gviz", "extdata", "gapped.bam"), isPaired = TRUE)
bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr1", start = afrom, end = ato, filter = list(with_ox_refseq_mrna = TRUE),stacking = "dense")
plotTracks(c(bmt, alTrack), from = afrom, to = ato,chromosome = "chr12")

aTrack <- AnnotationTrack(start = c(10, 40, 120), width = 15, chromosome = "chrX", strand = c("+", "*", "-"), id = c("Huey", "Dewey", "Louie"), genome = "ss10", name = "foo")
chanplotTracks(aTrack)
