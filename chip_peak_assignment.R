## loading packages

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
library(clusterProfiler)

##############################
#mouse analysis
##############################
mouse_genes = read.delim("./mouse/mouse_RPgenes.tsv")

#replicate 1
peakAnno_rep1 <- annotatePeak('./mouse/mouse_replicate1.bed', tssRegion=c(-3000, 3000), TxDb=txdb, annoDb = "org.Mm.eg.db")
all_peaks_mouse1 <- as.GRanges(peakAnno_rep1) #use GRanges rather than data.frame so that after filtering you can use covplot
filtered_peaks_mouse1 <- all_peaks_mouse1[all_peaks_mouse1$geneId %in% mouse_genes$Gene.ID, ]
head(filtered_peaks_mouse1)
unique(filtered_peaks_mouse1$SYMBOL)
filename <- file.path(getwd(), "./mouse/assigned_peaks/RPgenePeaks_mouse_rep1.csv")
write.csv(x=filtered_peaks_mouse1, file=filename, row.names = FALSE)

#replicate 2
peakAnno_rep2 <- annotatePeak('./mouse/mouse_replicate2.bed', tssRegion=c(-3000, 3000), TxDb=txdb, annoDb = "org.Mm.eg.db")
all_peaks_mouse2 <- as.GRanges(peakAnno_rep2) #use GRanges rather than data.frame so that after filtering you can use covplot
filtered_peaks_mouse2 <- all_peaks_mouse2[all_peaks_mouse2$geneId %in% mouse_genes$Gene.ID, ]
head(filtered_peaks_mouse2)
unique(filtered_peaks_mouse2$SYMBOL)
filename <- file.path(getwd(), "./mouse/assigned_peaks/RPgenePeaks_mouse_rep2.csv")
write.csv(x=filtered_peaks_mouse2, file=filename, row.names = FALSE)

plotAnnoPie(peakAnno_rep1)
dev.copy(png,'./mouse/graphs/annopie_mouse1.png')
dev.off()
plotAnnoPie(peakAnno_rep2)
dev.copy(png,'./mouse/graphs/annopie_mouse2.png')
dev.off()

plotAnnoBar(peakAnno_rep1)
dev.copy(png,'./mouse/graphs/annobar_mouse1.png', width=480, height=240)
dev.off()
plotAnnoBar(peakAnno_rep2)
dev.copy(png,'./mouse/graphs/annobar_mouse2.png', width=480, height=240)
dev.off()

plotDistToTSS(peakAnno_rep1)
dev.copy(png,'./mouse/graphs/DistToTSS_mouse1.png', width=480, height=240)
dev.off()
plotDistToTSS(peakAnno_rep2)
dev.copy(png,'./mouse/graphs/DistToTSS_mouse2.png', width=480, height=240)
dev.off()


replicate1_peak <- readPeakFile('./mouse/mouse_replicate1.bed')
replicate2_peak <- readPeakFile('./mouse/mouse_replicate2.bed')

covplot(replicate1_peak, weightCol = 'V5')
dev.copy(png,'./mouse/graphs/peak_coverage_mouse1.png')
dev.off()
covplot(replicate2_peak, weightCol = 'V5')
dev.copy(png,'./mouse/graphs/peak_coverage_mouse2.png')
dev.off()
covplot(filtered_peaks_mouse1, weightCol = 'V5')
dev.copy(png,'./mouse/graphs/peak_coverage_mouse1_RPgenes.png')
dev.off()
covplot(filtered_peaks_mouse2, weightCol = 'V5')
dev.copy(png,'./mouse/graphs/peak_coverage_mouse2_RPgenes.png')
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_replicate1 <- getTagMatrix(replicate1_peak, windows=promoter)
tagMatrix_replicate2 <- getTagMatrix(replicate2_peak, windows=promoter)
#tagHeatmap(tagMatrix_replicate1, xlim=c(-3000, 3000), color="red")
#dev.copy(png,'./mouse/graphs/heatmap_mouse1.png')
#dev.off()
#tagHeatmap(tagMatrix_replicate2, xlim=c(-3000, 3000), color="red")
#dev.copy(png,'./mouse/graphs/heatmap_mouse2.png')
#dev.off()

## get binding profile
plotAvgProf(tagMatrix_replicate1, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./mouse/graphs/avgprof_mouse1.png')
dev.off()
plotAvgProf(tagMatrix_replicate2, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./mouse/graphs/avgprof_mouse2.png')
dev.off()


