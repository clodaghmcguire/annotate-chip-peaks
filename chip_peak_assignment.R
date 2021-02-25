## loading packages

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
library(clusterProfiler)

## read in list of ribosomal protein genes
genes = read.delim("ribsome_geneset.txt")

## datasets used for PLOS paper
## datasets can be downloaded from GSE36994
## datasets are in hg18 but could be converted to hg19 or hg38 using UCSC liftover
## annotate peaks
peakAnno_GSM970258_A <- annotatePeak('GSM970258_GATA1-A_peaks.bed.gz', tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_GSM970257_F <- annotatePeak('GSM970257_GATA1-F_peaks.bed.gz', tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Hs.eg.db")
## get list of assigned peaks
all_peaks_GSM970258_A <- as.GRanges(peakAnno_GSM970258_A) #use GRanges rather than data.frame so that after filtering you can use covplot
head(all_peaks_GSM970258_A)
all_peaks_GSM970257_F <- as.GRanges(peakAnno_GSM970257_F)
head(all_peaks_GSM970257_F)

## filter annoPeaks to just ribosomal protein genes
filtered_peaks_GSM970258_A <- all_peaks_GSM970258_A[all_peaks_GSM970258_A$geneId %in% genes$gene_id, ]
head(filtered_peaks_GSM970258_A)
filtered_peaks_GSM970257_F <- all_peaks_GSM970257_F[all_peaks_GSM970257_F$geneId %in% genes$gene_id, ]
filename <- file.path(getwd(), "RPgenePeaks_GSM970258_A")
write.csv(x=filtered_peaks_GSM970258_A, file=filename)
filename <- file.path(getwd(), "RPgenePeaks_GSM970257_F")
write.csv(x=filtered_peaks_GSM970257_F, file=filename)

GSM970258_A_peak <- readPeakFile('GSM970258_GATA1-A_peaks.bed.gz')
GSM970257_F_peak <- readPeakFile('GSM970257_GATA1-F_peaks.bed.gz')

covplot(GSM970258_A_peak, weightCol = 'V5')
dev.copy(png,'peak_coverage_GSM970258_A.png')
dev.off()
covplot(filtered_peaks_GSM970258_A, weightCol = 'V5')
dev.copy(png,'peak_coverage_GSM970258_A_RPgenes.png')
dev.off()
covplot(GSM970257_F_peak, weightCol = 'V5')
dev.copy(png,'peak_coverage_GSM970257_F.png')
dev.off()
covplot(filtered_peaks_GSM970257_F, weightCol = 'V5')
dev.copy(png,'peak_coverage_GSM970257_F_RPgenes.png')
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=10000, downstream=10000)
tagMatrix_GSM970258_A <- getTagMatrix(GSM970258_A_peak, windows=promoter)
tagMatrix_GSM970257_F <- getTagMatrix(GSM970257_F_peak, windows=promoter)
tagHeatmap(tagMatrix_GSM970258_A, xlim=c(-10000, 10000), color="red")
dev.copy(png,'heatmap_GSM970258_A.png')
dev.off()
tagHeatmap(tagMatrix_GSM970257_F, xlim=c(-10000, 10000), color="red")
dev.copy(png,'heatmap_GSM970257_F.png')
dev.off()

tagMatrix_RP_A <- getTagMatrix(filtered_peaks_GSM970258_A, windows=promoter)
tagMatrix_RP_F <- getTagMatrix(filtered_peaks_GSM970257_F, windows=promoter)

## get binding profile
plotAvgProf(tagMatrix_GSM970258_A, xlim=c(-10000, 10000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'avgprof_GSM970258_A.png')
dev.off()
plotAvgProf(tagMatrix_GSM970257_F, xlim=c(-10000, 10000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'avgprof_GSM970257_F.png')
dev.off()
plotAvgProf(tagMatrix_RP_A, xlim=c(-10000, 10000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'avgprof_RP_A.png')
dev.off()
plotAvgProf(tagMatrix_RP_F, xlim=c(-10000, 10000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'avgprof_RP_F.png')
dev.off()

#peak annotation
plotAnnoPie(peakAnno_GSM970258_A)
dev.copy(png,'annopie_GSM970258_A.png')
dev.off()
plotAnnoPie(peakAnno_GSM970257_F)
dev.copy(png,'annopie_GSM970257_F.png')
dev.off()
plotAnnoPie(peakAnno_RP_A)
dev.copy(png,'annopie_RP_A.png')
dev.off()
plotAnnoPie(peakAnno_RP_F)
dev.copy(png,'annopie_RP_F.png')
dev.off()


#TF binding
plotAnnoBar(peakAnno_GSM970258_A)
dev.copy(png,'annobar_GSM970258_A.png')
dev.off()
plotAnnoBar(peakAnno_GSM970257_F)
dev.copy(png,'annobar_GSM970257_F.png')
dev.off()
plotAnnoBar(peakAnno_RP_A)
dev.copy(png,'annobar_RP_A.png')
dev.off()
plotAnnoBar(peakAnno_RP_F)
dev.copy(png,'annobar_RP_F.png')
dev.off()


#distribution of TF binding relative to TSS
plotDistToTSS(peakAnno_GSM970258_A)
dev.copy(png,'DistToTSS_GSM970258_A.png')
dev.off()
plotDistToTSS(peakAnno_GSM970257_F)
dev.copy(png,'DistToTSS_GSM970257_F.png')
dev.off()
plotDistToTSS(peakAnno_RP_A)
dev.copy(png,'DistToTSS_RP_A.png')
dev.off()
plotDistToTSS(peakAnno_RP_F)
dev.copy(png,'DistToTSS_RP_F.png')
dev.off()