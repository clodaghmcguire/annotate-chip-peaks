library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

## read in list of ribosomal protein genes
genes = read.delim("ribsome_geneset.txt")

## datasets can be downloaded from GSE36994
## original FASTQ aligned to hg19 and peaks called using macs3
## annotate peaks
peakAnno_GATA1 <- annotatePeak('./GSE36994_hg19_macs3/GATA1_macs3_peaks_enrichment5.bed', tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Hs.eg.db")

## get list of assigned peaks
all_peaks_GATA1 <- as.GRanges(peakAnno_GATA1) #use GRanges rather than data.frame so that after filtering you can use covplot
head(all_peaks_GATA1)

## filter annoPeaks to just ribosomal protein genes
filtered_peaks_GATA1 <- all_peaks_GATA1[all_peaks_GATA1$geneId %in% genes$gene_id, ]
head(filtered_peaks_GATA1)

filename <- file.path(getwd(), "./GSE36994_hg19_macs3/assigned_peaks/GATA1_chipseeker.csv")
write.csv(x=filtered_peaks_GATA1, file=filename, row.names = FALSE)


GATA1_peak <- readPeakFile('./GSE36994_hg19_macs3/GATA1_macs3_peaks_enrichment5.bed')


covplot(GATA1_peak)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/peak_coverage_GATA1.png')
dev.off()
covplot(filtered_peaks_GATA1)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/peak_coverage_GATA1_RPgenes.png')
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_GATA1 <- getTagMatrix(GATA1_peak, windows=promoter)


## get binding profile
plotAvgProf(tagMatrix_GATA1, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./GSE36994_hg19_macs3/graphs/avgprof_GATA1.png')
dev.off()

#peak annotation
plotAnnoPie(peakAnno_GATA1)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/annopie_GATA1.png')
dev.off()
plotAnnoPie(filtered_peaks_GATA1)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/annopie_GATA1_RP.png')
dev.off()


#TF binding
plotAnnoBar(peakAnno_GATA1)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/annobar_GATA1.png', width=480, height=240)
dev.off()


#distribution of TF binding relative to TSS
plotDistToTSS(peakAnno_GATA1, title = "Distribution of TF binding relative to TSS in adult cells")
dev.copy(png,'./GSE36994_hg19_macs3/graphs/DistToTSS_GATA1.png', width=480, height=240)
dev.off()


#### comparison against original bed file (macs and score >1000) ####
library(VennDiagram)
library(RColorBrewer)
GSE36994_original <- read.csv('./GSE36994_hg19/assigned_peaks/RPgenePeaks_GSM970258_A_hg19_chipseeker.csv')

venn.diagram(x = list(unique(GSE36994_original$SYMBOL), unique(filtered_peaks_GATA1$SYMBOL)), category.names = c("Original", "Re-analysis"), filename = './GSE36994_hg19_macs3/dataset_comparison/original_vs_macs3_re-analysis.png',
             output=TRUE, imagetype="png"
)

#### comparison of all datasets: ENCODE, GSE36994 and GSE43625 ####
ENCODE_GATA1 <- read.csv('./ENCODE/assigned_peaks/RPgenePeaks_ENCODE_A_hg19.csv')
GSE43625 <- read.csv('./GSE43625_hg19/assigned_peaks/GATA1.csv')

venn.diagram(x = list(unique(GSE36994_original$SYMBOL), unique(ENCODE_GATA1$SYMBOL), unique(GSE43625$SYMBOL)), category.names = c("GSE36994", "ENCODE", "GSE43625"), filename = './GSE36994_hg19_macs3/dataset_comparison/three_datasets_overlap.png',
             output=TRUE, imagetype="png"
)

venn.diagram(x = list(unique(GSE36994_original$SYMBOL), unique(ENCODE_GATA1$SYMBOL), unique(GSE43625$SYMBOL), unique(filtered_peaks_GATA1$SYMBOL)), category.names = c("GSE36994", "ENCODE", "GSE43625", "GSE36994 re-analysis"), filename = './GSE36994_hg19_macs3/dataset_comparison/all_datasets_overlap.png',
             output=TRUE, imagetype="png"
)

unique(ENCODE_GATA1$SYMBOL)
unique(GSE36994_original$SYMBOL)
unique(GSE43625$SYMBOL)
unique(filtered_peaks_GATA1$SYMBOL)
RPgenelist <- read.delim('RPgenes', header = TRUE, sep = '\t')
setdiff(unique(RPgenelist$Symbol), unique(filtered_peaks_GATA1$SYMBOL))
