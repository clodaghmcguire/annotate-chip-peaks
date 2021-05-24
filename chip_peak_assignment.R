## loading packages

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

## read in list of ribosomal protein genes
genes = read.delim("ribsome_geneset.txt")

## datasets used for PLOS paper
## datasets can be downloaded from GSE36994
## original datasets are in hg18 but were converted to hg19 UCSC liftover
## annotate peaks
peakAnno_GSM970258_A <- annotatePeak('./GSE36994_hg19/GSM970258_GATA1-A_hg19.bed', tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_GSM970257_F <- annotatePeak('./GSE36994_hg19/GSM970257_GATA1-F_hg19.bed', tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Hs.eg.db")
## get list of assigned peaks
all_peaks_GSM970258_A <- as.GRanges(peakAnno_GSM970258_A) #use GRanges rather than data.frame so that after filtering you can use covplot
head(all_peaks_GSM970258_A)
all_peaks_GSM970257_F <- as.GRanges(peakAnno_GSM970257_F)
head(all_peaks_GSM970257_F)

## filter annoPeaks to just ribosomal protein genes
filtered_peaks_GSM970258_A <- all_peaks_GSM970258_A[all_peaks_GSM970258_A$geneId %in% genes$gene_id, ]
head(filtered_peaks_GSM970258_A)
filtered_peaks_GSM970257_F <- all_peaks_GSM970257_F[all_peaks_GSM970257_F$geneId %in% genes$gene_id, ]
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/RPgenePeaks_GSM970258_A_hg19_chipseeker.csv")
write.csv(x=filtered_peaks_GSM970258_A, file=filename, row.names = FALSE)
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/RPgenePeaks_GSM970257_F_hg19.csv")
write.csv(x=filtered_peaks_GSM970257_F, file=filename, row.names = FALSE)

GSM970258_A_peak <- readPeakFile('./GSE36994_hg19/GSM970258_GATA1-A_hg19.bed')
GSM970257_F_peak <- readPeakFile('./GSE36994_hg19/GSM970257_GATA1-F_hg19.bed')

covplot(GSM970258_A_peak, weightCol = 'V5')
dev.copy(png,'./GSE36994_hg19/graphs/peak_coverage_GSM970258_A.png')
dev.off()
covplot(filtered_peaks_GSM970258_A, weightCol = 'V5')
dev.copy(png,'./GSE36994_hg19/graphs/peak_coverage_GSM970258_A_RPgenes.png')
dev.off()
covplot(GSM970257_F_peak, weightCol = 'V5')
dev.copy(png,'./GSE36994_hg19/graphs/peak_coverage_GSM970257_F.png')
dev.off()
covplot(filtered_peaks_GSM970257_F, weightCol = 'V5')
dev.copy(png,'./GSE36994_hg19/graphs/peak_coverage_GSM970257_F_RPgenes.png')
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_GSM970258_A <- getTagMatrix(GSM970258_A_peak, windows=promoter)
tagMatrix_GSM970257_F <- getTagMatrix(GSM970257_F_peak, windows=promoter)
#tagHeatmap(tagMatrix_GSM970258_A, xlim=c(-10000, 10000), color="red")
#dev.copy(png,'./GSE36994_hg19/graphs/heatmap_GSM970258_A.png')
#dev.off()
#tagHeatmap(tagMatrix_GSM970257_F, xlim=c(-10000, 10000), color="red")
#dev.copy(png,'./GSE36994_hg19/graphs/heatmap_GSM970257_F.png')
#dev.off()

## get binding profile
plotAvgProf(tagMatrix_GSM970258_A, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./GSE36994_hg19/graphs/avgprof_GSM970258_A.png')
dev.off()
plotAvgProf(tagMatrix_GSM970257_F, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./GSE36994_hg19/graphs/avgprof_GSM970257_F.png')
dev.off()

#peak annotation
plotAnnoPie(peakAnno_GSM970258_A)
dev.copy(png,'./GSE36994_hg19/graphs/annopie_GSM970258_A.png')
dev.off()
plotAnnoPie(peakAnno_GSM970257_F)
dev.copy(png,'./GSE36994_hg19/graphs/annopie_GSM970257_F.png')
dev.off()


#TF binding
plotAnnoBar(peakAnno_GSM970258_A)
dev.copy(png,'./GSE36994_hg19/graphs/annobar_GSM970258_A.png', width=480, height=240)
dev.off()
plotAnnoBar(peakAnno_GSM970257_F)
dev.copy(png,'./GSE36994_hg19/graphs/annobar_GSM970257_F.png', width=480, height=240)
dev.off()


#distribution of TF binding relative to TSS
plotDistToTSS(peakAnno_GSM970258_A, title = "Distribution of TF binding relative to TSS in adult cells")
dev.copy(png,'./GSE36994_hg19/graphs/DistToTSS_GSM970258_A.png', width=480, height=240)
dev.off()
plotDistToTSS(peakAnno_GSM970257_F, title = "Distribution of TF binding relative to TSS in fetal cells")
dev.copy(png,'./GSE36994_hg19/graphs/DistToTSS_GSM970257_F.png', width=480, height=240)
dev.off()

##comparison of adult vs fetal
files <- list(Adult = './GSE36994_hg19/GSM970258_GATA1-A_hg19.bed', Fetal = './GSE36994_hg19/GSM970257_GATA1-F_hg19.bed')
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./GSE36994_hg19/graphs/avgprof_GATA1.png')
dev.off()
PlotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95, resample=500, facet="row")
dev.copy(png,'./GSE36994_hg19/graphs/avgprof_conf_GATA1.png')
dev.off()
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
dev.copy(png,'./GSE36994_hg19/graphs/annobar_GATA1.png', width=480, height=240)
dev.off()
plotDistToTss(peakAnnoList)
dev.copy(png,'./GSE36994_hg19/graphs/DistToTSS_GATA1.png', width=480, height=240)
dev.off()

############################
## chromotin marks ##
adult_files <- list(GATA1 = './GSE36994_hg19/GSM970258_GATA1-A_hg19.bed', 
                    TAL1 = './GSE36994_hg19/GSM908055_TAL1-A_hg19.bed',
                    PolII = './GSE36994_hg19/GSM908069_PolII-A_peaks_hg19.bed',
                    H3K27ac = './GSE36994_hg19/GSM908051_H3K27ac-A_peaks_hg19.bed',
                    H3K36me3 = './GSE36994_hg19/GSM908047_H3K36me3-A_peaks_hg19.bed',
                    H3K27me3 = './GSE36994_hg19/GSM908043_H3K27me3-A_peaks_hg19.bed',
                    H3K9me3 = './GSE36994_hg19/GSM908041_H3K9me3-A_peaks_hg19.bed',
                    H3K4me3 = './GSE36994_hg19/GSM908039_H3K4me3-A_peaks_hg19.bed',
                    H3K4me1 = './GSE36994_hg19/GSM908035_H3K4me1-A_peaks_hg19.bed')

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

tagMatrixList <- lapply(adult_files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), facet="row")
par(mar=c(1,1,1,1))
#tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
peakAnnoList <- lapply(adult_files, annotatePeak, TxDb=txdb, annoDb="org.Hs.eg.db",
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

#all_genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
rp_gene_peaks <- lapply(peakAnnoList, function(i) as.data.frame(i)[as.data.frame(i)$geneId %in% genes$gene_id, ])
sapply(names(rp_gene_peaks), function (x) write.csv(rp_gene_peaks[[x]], file=paste(paste("./GSE36994_hg19/assigned_peaks/histones/", x), "csv", sep="."), row.names = FALSE )   )
rp_genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId[as.data.frame(i)$geneId %in% genes$gene_id])
vennplot(rp_genes)
dev.copy(png,'./GSE36994_hg19/graphs/histone_gene_overlap.png')
dev.off()

#### TAL1 ####
## datasets used for PLOS paper
## datasets can be downloaded from GSE36994
## datasets are in hg18 but were converted to hg19 or hg38 using UCSC liftover
## annotate peaks
peakAnno_TAL1_A <- annotatePeak('./GSE36994_hg19/GSM908055_TAL1-A_hg19.bed', tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_TAL1_F <- annotatePeak('./GSE36994_hg19/GSM908054_TAL1-F_hg19.bed', tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
## get list of assigned peaks
all_peaks_TAL1_A <- as.GRanges(peakAnno_TAL1_A) #use GRanges rather than data.frame so that after filtering you can use covplot
head(all_peaks_TAL1_A)
all_peaks_TAL1_F <- as.GRanges(peakAnno_TAL1_F)
head(all_peaks_TAL1_F)

## filter annoPeaks to just ribosomal protein genes
filtered_peaks_TAL1_A <- all_peaks_TAL1_A[all_peaks_TAL1_A$geneId %in% genes$gene_id, ]
head(filtered_peaks_TAL1_A)
filtered_peaks_TAL1_F <- all_peaks_TAL1_F[all_peaks_TAL1_F$geneId %in% genes$gene_id, ]
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/RPgenePeaks_TAL1_A_hg19.csv")
write.csv(x=filtered_peaks_TAL1_A, file=filename, row.names = FALSE)
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/RPgenePeaks_TAL1_F_hg19.csv")
write.csv(x=filtered_peaks_TAL1_F, file=filename, row.names = FALSE)

TAL1_A_peak <- readPeakFile('./GSE36995_hg19/TAL1-A_hg19.bed')
TAL1_F_peak <- readPeakFile('./GSE36995_hg19/TAL1-F_hg19.bed')

covplot(TAL1_A_peak)
dev.copy(png,'./GSE36994_hg19/graphs/peak_coverage_TAL1_A.png')
dev.off()
covplot(filtered_peaks_TAL1_A)
dev.copy(png,'./GSE36994_hg19/graphs/peak_coverage_TAL1_A_RPgenes.png')
dev.off()
covplot(TAL1_F_peak, weightCol = 'V5')
dev.copy(png,'./GSE36994_hg19/graphs/peak_coverage_TAL1_F.png')
dev.off()
covplot(filtered_peaks_TAL1_F, weightCol = 'V5')
dev.copy(png,'./GSE36994_hg19/graphs/peak_coverage_TAL1_F_RPgenes.png')
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_TAL1_A <- getTagMatrix(TAL1_A_peak, windows=promoter)
tagMatrix_TAL1_F <- getTagMatrix(TAL1_F_peak, windows=promoter)
#tagHeatmap(tagMatrix_TAL1_A, xlim=c(-3000, 3000), color="red")
#dev.copy(png,'./GSE36994_hg19/heatmap_TAL1_A.png')
#dev.off()
#tagHeatmap(tagMatrix_GSM970257_F, xlim=c(-3000, 3000), color="red")
#dev.copy(png,'./GSE36994_hg19/heatmap_TAL1_F.png')
#dev.off()

## get binding profile
plotAvgProf(tagMatrix_TAL1_A, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./GSE36994_hg19/graphs/avgprof_TAL18_A.png')
dev.off()
plotAvgProf(tagMatrix_TAL1_F, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./GSE36994_hg19/graphs/avgprof_TAL1_F.png')
dev.off()

#peak annotation
plotAnnoPie(peakAnno_TAL1_A)
dev.copy(png,'./GSE36994_hg19/graphs/annopie_TAL1_A.png')
dev.off()
plotAnnoPie(peakAnno_TAL1_F)
dev.copy(png,'./GSE36994_hg19/graphs/annopie_TAL1_F.png')
dev.off()

#TF binding
plotAnnoBar(peakAnno_TAL1_A)
dev.copy(png,'./GSE36994_hg19/graphs/annobar_TAL1_A.png', width=600, height=300)
dev.off()
plotAnnoBar(peakAnno_TAL1_F)
dev.copy(png,'./GSE36994_hg19/graphs/annobar_TAL1_F.png', width=600, height=300)
dev.off()


#distribution of TF binding relative to TSS
plotDistToTSS(peakAnno_TAL1_A)
dev.copy(png,'./GSE36994_hg19/graphs/DistToTSS_TAL1_A.png', width=600, height=300)
dev.off()
plotDistToTSS(peakAnno_TAL1_F)
dev.copy(png,'./GSE36994_hg19/graphs/DistToTSS_TAL1_F.png', width=600, height=300)
dev.off()


##########################################################################
##### Peak assignment for GSE43625 KLF1 #####
##########################################################################
files <- list(GATA1 = './GSE43625_hg19/GSM1067274_Erythroid_GATA1_peaks_hg19.bed.gz', 
              TAL1 = './GSE43625_hg19/GSM1067277_Erythroid_TAL1_peaks_hg19.bed.gz',
              KLF1 = './GSE43625_hg19/GSM1067275_Erythroid_KLF1_peaks_hg19.bed.gz',
              NFE2 = './GSE43625_hg19/GSM1067276_Erythroid_NFE2_peaks_hg19.bed.gz'
)

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

tagMatrixList_files <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList_files, xlim=c(-3000, 3000))
dev.copy(png,'./GSE43625_hg19/graphs/avgprof.png', width=600, height=300)
dev.off()
plotAvgProf(tagMatrixList_files, xlim=c(-3000, 3000), facet="row")
dev.copy(png,'./GSE43625_hg19/graphs/avgprof_row.png', width=600, height=300)
dev.off()

peakAnnoList_files <- lapply(files, annotatePeak, TxDb=txdb, annoDb="org.Hs.eg.db",
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList_files)
dev.copy(png,'./GSE43625_hg19/graphs/annobar.png', width=600, height=300)
dev.off()
plotDistToTSS(peakAnnoList_files)
dev.copy(png,'./GSE43625_hg19/graphs/DistToTss.png', width=600, height=300)
dev.off()

rp_gene_peaks <- lapply(peakAnnoList_files, function(i) as.data.frame(i)[as.data.frame(i)$geneId %in% genes$gene_id, ])
sapply(names(rp_gene_peaks), function (x) write.csv(rp_gene_peaks[[x]], file=paste(paste("./GSE43625_hg19/assigned_peaks/", x), "csv", sep="."), row.names = FALSE )   )
sapply(names(rp_gene_peaks), function (x) write.table(rp_gene_peaks[[x]], file=paste(paste("./GSE43625_hg19/assigned_peaks/", x), "bed", sep="."), sep = "\t", row.names = FALSE )   )

rp_genes <- lapply(rp_gene_peaks, function(i) unique(as.data.frame(i)$geneId))
vennplot(rp_genes)
dev.copy(png,'./GSE43625_hg19/graphs/gene_overlap.png')
dev.off()

###############################################################################
#### ENCODE ####
###############################################################################
peakAnno_ENCODE_adult <- annotatePeak('./ENCODE/ENCFF957CWW.bed.gz', tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_ENCODE_fetal <- annotatePeak('./ENCODE/ENCFF762DRA.bed.gz', tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

## get list of assigned peaks
all_peaks_ENCODE_adult <- as.GRanges(peakAnno_ENCODE_adult) #use GRanges rather than data.frame so that after filtering you can use covplot
all_peaks_ENCODE_fetal <- as.GRanges(peakAnno_ENCODE_fetal)

## filter annoPeaks to just ribosomal protein genes
filtered_peaks_ENCODE_adult <- all_peaks_ENCODE_adult[all_peaks_ENCODE_adult$geneId %in% genes$gene_id, ]
filename <- file.path(getwd(), "./ENCODE/assigned_peaks/RPgenePeaks_ENCODE_A_hg19.csv")
write.csv(x=filtered_peaks_ENCODE_adult, file=filename, row.names = FALSE)

filtered_peaks_ENCODE_fetal <- all_peaks_ENCODE_fetal[all_peaks_ENCODE_fetal$geneId %in% genes$gene_id, ]
filename <- file.path(getwd(), "./ENCODE/assigned_peaks/RPgenePeaks_ENCODE_F_hg19.csv")
write.csv(x=filtered_peaks_ENCODE_fetal, file=filename, row.names = FALSE)

ENCODE_peak_adult <- readPeakFile('./ENCODE/ENCFF957CWW.bed.gz')
ENCODE_peak_fetal <- readPeakFile('./ENCODE/ENCFF762DRA.bed.gz')

covplot(ENCODE_peak_adult, weightCol = 'V5')
dev.copy(png,'./ENCODE/graphs/peak_coverage_ENCODE_A.png')
dev.off()
covplot(filtered_peaks_ENCODE_adult, weightCol = 'V5')
dev.copy(png,'./ENCODE/graphs/peak_coverage_ENCODE_RPgenes_A.png')
dev.off()

covplot(ENCODE_peak_fetal, weightCol = 'V5')
dev.copy(png,'./ENCODE/graphs/peak_coverage_ENCODE_F.png')
dev.off()
covplot(filtered_peaks_ENCODE_fetal, weightCol = 'V5')
dev.copy(png,'peak_coverage_ENCODE_RPgenes_F.png')
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_ENCODE_adult <- getTagMatrix(ENCODE_peak_adult, windows=promoter)
tagMatrix_ENCODE_fetal <- getTagMatrix(ENCODE_peak_fetal, windows=promoter)

## get binding profile
plotAvgProf(tagMatrix_ENCODE_adult, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./ENCODE/graphs/avgprof_ENCODE_A.png')
dev.off()
plotAvgProf(tagMatrix_ENCODE_fetal, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./ENCODE/graphs/avgprof_ENCODE_F.png')
dev.off()

#peak annotation
plotAnnoPie(peakAnno_ENCODE_adult)
dev.copy(png,'./ENCODE/graphs/annopie_ENCODE_A.png')
dev.off()
plotAnnoPie(peakAnno_ENCODE_fetal)
dev.copy(png,'./ENCODE/graphs/annopie_ENCODE_F.png')
dev.off()

#TF binding
plotAnnoBar(peakAnno_ENCODE_adult)
dev.copy(png,'./ENCODE/graphs/annobar_ENCODE_A.png', width=480, height=240)
dev.off()
plotAnnoBar(peakAnno_ENCODE_fetal)
dev.copy(png,'./ENCODE/graphs/annobar_ENCODE_F.png', width=480, height=240)
dev.off()

#distribution of TF binding relative to TSS
plotDistToTSS(peakAnno_ENCODE_adult, title = "Distribution of TF binding relative to TSS in adult cells")
dev.copy(png,'./ENCODE/graphs/DistToTSS_ENCODE_A.png', width=480, height=240)
dev.off()
plotDistToTSS(peakAnno_ENCODE_fetal, title = "Distribution of TF binding relative to TSS in adult cells")
dev.copy(png,'./ENCODE/graphs/DistToTSS_ENCODE_F.png', width=480, height=240)
dev.off()


##############################
#mouse analysis
##############################
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

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


