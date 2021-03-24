## loading packages

library(ChIPseeker)
#library(TxDb.Hsapiens.UCSC.hg18.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

## read in list of ribosomal protein genes
genes = read.delim("ribsome_geneset.txt")

## datasets used for PLOS paper
## datasets can be downloaded from GSE36994
## datasets are in hg18 but could be converted to hg19 or hg38 using UCSC liftover
## annotate peaks
peakAnno_GSM970258_A <- annotatePeak('GSM970258_GATA1-A_hg19.bed', tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_GSM970257_F <- annotatePeak('GSM970257_GATA1-F_hg19.bed', tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Hs.eg.db")
## get list of assigned peaks
all_peaks_GSM970258_A <- as.GRanges(peakAnno_GSM970258_A) #use GRanges rather than data.frame so that after filtering you can use covplot
head(all_peaks_GSM970258_A)
all_peaks_GSM970257_F <- as.GRanges(peakAnno_GSM970257_F)
head(all_peaks_GSM970257_F)

## filter annoPeaks to just ribosomal protein genes
filtered_peaks_GSM970258_A <- all_peaks_GSM970258_A[all_peaks_GSM970258_A$geneId %in% genes$gene_id, ]
head(filtered_peaks_GSM970258_A)
filtered_peaks_GSM970257_F <- all_peaks_GSM970257_F[all_peaks_GSM970257_F$geneId %in% genes$gene_id, ]
#filename <- file.path(getwd(), "RPgenePeaks_GSM970258_A")
filename <- file.path(getwd(), "RPgenePeaks_GSM970258_A_hg19")
write.csv(x=filtered_peaks_GSM970258_A, file=filename)
#filename <- file.path(getwd(), "RPgenePeaks_GSM970257_F")
filename <- file.path(getwd(), "RPgenePeaks_GSM970257_F_hg19")
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

##comparison of adult vs fetal
files <- list(Adult = 'GSM970258_GATA1-A_hg19.bed', Fetal = 'GSM970257_GATA1-F_hg19.bed')
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-10000, 10000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'avgprof.png')
dev.off()
PlotAvgProf(tagMatrixList, xlim=c(-10000, 10000), conf=0.95, resample=500, facet="row")
dev.copy(png,'avgprof_conf.png')
dev.off()
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-10000, 10000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
dev.copy(png,'annobar.png')
dev.off()
plotDistToTss(peakAnnoList)
dev.copy(png,'DistToTSS.png')
dev.off()


##############################
#mouse analysis
##############################
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
library(clusterProfiler)

mouse_genes = read.delim("mouse_RPgenes.tsv")

#replicate 1
peakAnno_rep1 <- annotatePeak('mouse_replicate1.bed', tssRegion=c(-10000, 10000), TxDb=txdb, annoDb = "org.Mm.eg.db")
all_peaks_mouse1 <- as.GRanges(peakAnno_rep1) #use GRanges rather than data.frame so that after filtering you can use covplot
filtered_peaks_mouse1 <- all_peaks_mouse1[all_peaks_mouse1$geneId %in% mouse_genes$Gene.ID, ]
head(filtered_peaks_mouse1)
unique(filtered_peaks_mouse1$SYMBOL)
filename <- file.path(getwd(), "RPgenePeaks_mouse_rep1")
write.csv(x=filtered_peaks_mouse1, file=filename)

#replicate 2
peakAnno_rep2 <- annotatePeak('mouse_replicate2.bed', tssRegion=c(-10000, 10000), TxDb=txdb, annoDb = "org.Mm.eg.db")
all_peaks_mouse2 <- as.GRanges(peakAnno_rep2) #use GRanges rather than data.frame so that after filtering you can use covplot
filtered_peaks_mouse2 <- all_peaks_mouse2[all_peaks_mouse2$geneId %in% mouse_genes$Gene.ID, ]
head(filtered_peaks_mouse2)
unique(filtered_peaks_mouse2$SYMBOL)
filename <- file.path(getwd(), "RPgenePeaks_mouse_rep2")
write.csv(x=filtered_peaks_mouse2, file=filename)

plotAnnoPie(peakAnno_rep1)
dev.copy(png,'annopie_mouse1.png')
dev.off()
plotAnnoPie(peakAnno_rep2)
dev.copy(png,'annopie_mouse2.png')
dev.off()

plotAnnoBar(peakAnno_rep1)
dev.copy(png,'annobar_mouse1.png')
dev.off()
plotAnnoBar(peakAnno_rep2)
dev.copy(png,'annobar_mouse2.png')
dev.off()

plotDistToTSS(peakAnno_rep1)
dev.copy(png,'DistToTSS_mouse1.png')
dev.off()
plotDistToTSS(peakAnno_rep2)
dev.copy(png,'DistToTSS_mouse2.png')
dev.off()


replicate1_peak <- readPeakFile('mouse_replicate1.bed')
replicate2_peak <- readPeakFile('mouse_replicate2.bed')

covplot(replicate1_peak, weightCol = 'V5')
dev.copy(png,'peak_coverage_mouse1.png')
dev.off()
covplot(replicate2_peak, weightCol = 'V5')
dev.copy(png,'peak_coverage_mouse2.png')
dev.off()
covplot(filtered_peaks_mouse1, weightCol = 'V5')
dev.copy(png,'peak_coverage_mouse1_RPgenes.png')
dev.off()
covplot(filtered_peaks_mouse2, weightCol = 'V5')
dev.copy(png,'peak_coverage_mouse2_RPgenes.png')
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=10000, downstream=10000)
tagMatrix_replicate1 <- getTagMatrix(replicate1_peak, windows=promoter)
tagMatrix_replicate2 <- getTagMatrix(replicate2_peak, windows=promoter)
tagHeatmap(tagMatrix_replicate1, xlim=c(-10000, 10000), color="red")
dev.copy(png,'heatmap_mouse1.png')
dev.off()
tagHeatmap(tagMatrix_replicate2, xlim=c(-10000, 10000), color="red")
dev.copy(png,'heatmap_mouse2.png')
dev.off()

## get binding profile
plotAvgProf(tagMatrix_replicate1, xlim=c(-10000, 10000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'avgprof_mouse1.png')
dev.off()
plotAvgProf(tagMatrix_replicate2, xlim=c(-10000, 10000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'avgprof_mouse2.png')
dev.off()



############################
## chromotin marks ##
adult_files <- list(GATA1 = 'GSM970258_GATA1-A_hg19.bed', 
                    PolII = 'GSM908069_PolII-A_peaks_hg19.bed',
                    H3K27ac = 'GSM908051_H3K27ac-A_peaks_hg19.bed',
                    #H3K9ac = 'GSM908049_H3K9ac-A_peaks_hg19.bed',
                    H3K36me3 = 'GSM908047_H3K36me3-A_peaks_hg19.bed',
                    #H3K36me2 = 'GSM908045_H3K36me2-A_peaks_hg19.bed',
                    H3K27me3 = 'GSM908043_H3K27me3-A_peaks_hg19.bed',
                    H3K9me3 = 'GSM908041_H3K9me3-A_peaks_hg19.bed',
                    H3K4me3 = 'GSM908039_H3K4me3-A_peaks_hg19.bed',
                    #H3K4me2 = 'GSM908037_H3K4me2-A_peaks_hg19.bed',
                    H3K4me1 = 'GSM908035_H3K4me1-A_peaks_hg19.bed')

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

tagMatrixList <- lapply(adult_files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), facet="row")
par(mar=c(1,1,1,1))
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
peakAnnoList <- lapply(adult_files, annotatePeak, TxDb=txdb, annoDb="org.Hs.eg.db",
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

par(mar=c(1,1,1,1))
for(i in length(tagMatrixList)) {
 tagHeatmap(
    tagMatrixList[[i]],
    xlim = c(-3000, 3000)
  )
}
dev.copy(png,'histone_test.png')
dev.off()
  #tagMatrix <- getTagMatrix(adult_files[[i]], windows=promoter)
  #plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
  #dev.copy(png,'avgprof_',chip_names[i],'.png')
  #dev.off()
  #tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color=NULL)
  #dev.copy(png, 'heatmap_',chip_names[i],'.png')
  #dev.off()
  #peakAnno <- annotatePeak(adult_files[[i]], tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
  #plotAnnoBar(peakAnno)
  #dev.copy(png,'annobar_',chip_names[i],'.png')
  #dev.off()
  #plotDistToTss(peakAnno)
  #dev.copy(png,'DistToTSS_',chip_names[i],'.png')
  #dev.off()
  
  #all_peaks <- as.GRanges(peakAnno)
  #filtered_peaks <- all_peaks[all_peaks$geneId %in% genes$gene_id, ]
  #filename <- file.path(getwd(), chip_names[i])
  #write.csv(x=filtered_peaks, file=filename)


#all_genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
rp_gene_peaks <- lapply(peakAnnoList, function(i) as.data.frame(i)[as.data.frame(i)$geneId %in% genes$gene_id, ])
sapply(names(rp_gene_peaks), function (x) write.table(rp_gene_peaks[[x]], file=paste(x, "txt", sep=".") )   )
#rp_genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId[as.data.frame(i)$geneId %in% genes$gene_id])
#vennplot(rp_genes)

rp_gene_peaks <- lapply(peakAnnoList, function(i) as.GRanges(i)[as.GRanges(i)$geneId %in% genes$gene_id, ])
tagMatrixListFiltered <- lapply(rp_gene_peaks, getTagMatrix, windows=promoter)
par(mar=c(1,1,1,1))
tagHeatmap(tagMatrixListFiltered, xlim=c(-3000, 3000), color=NULL)
dev.copy(png,'H3K27ac.png')
dev.off()


#### TAL1 ####
## datasets used for PLOS paper
## datasets can be downloaded from GSE36994
## datasets are in hg18 but could be converted to hg19 or hg38 using UCSC liftover
## annotate peaks
peakAnno_TAL1_A <- annotatePeak('TAL1-A_hg19.bed', tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_TAL1_F <- annotatePeak('TAL1-F_hg19.bed', tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Hs.eg.db")
## get list of assigned peaks
all_peaks_TAL1_A <- as.GRanges(peakAnno_TAL1_A) #use GRanges rather than data.frame so that after filtering you can use covplot
head(all_peaks_TAL1_A)
all_peaks_TAL1_F <- as.GRanges(peakAnno_TAL1_F)
head(all_peaks_TAL1_F)

## filter annoPeaks to just ribosomal protein genes
filtered_peaks_TAL1_A <- all_peaks_TAL1_A[all_peaks_TAL1_A$geneId %in% genes$gene_id, ]
head(filtered_peaks_TAL1_A)
filtered_peaks_TAL1_F <- all_peaks_TAL1_F[all_peaks_TAL1_F$geneId %in% genes$gene_id, ]
filename <- file.path(getwd(), "RPgenePeaks_TAL1_A_hg19")
write.csv(x=filtered_peaks_TAL1_A, file=filename)
filename <- file.path(getwd(), "RPgenePeaks_TAL1_F_hg19")
write.csv(x=filtered_peaks_TAL1_F, file=filename)

TAL1_A_peak <- readPeakFile('TAL1-A_hg19.bed')
TAL1_F_peak <- readPeakFile('TAL1-F_hg19.bed')

covplot(TAL1_A_peak, weightCol = 'V5')
dev.copy(png,'peak_coverage_TAL1_A.png')
dev.off()
covplot(filtered_peaks_TAL1_A, weightCol = 'V5')
dev.copy(png,'peak_coverage_TAL1_A_RPgenes.png')
dev.off()
covplot(TAL1_F_peak, weightCol = 'V5')
dev.copy(png,'peak_coverage_TAL1_F.png')
dev.off()
covplot(filtered_peaks_TAL1_F, weightCol = 'V5')
dev.copy(png,'peak_coverage_TAL1_F_RPgenes.png')
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_TAL1_A <- getTagMatrix(TAL1_A_peak, windows=promoter)
tagMatrix_TAL1_F <- getTagMatrix(TAL1_F_peak, windows=promoter)
tagHeatmap(tagMatrix_TAL1_A, xlim=c(-3000, 3000), color="red")
dev.copy(png,'heatmap_TAL1_A.png')
dev.off()
tagHeatmap(tagMatrix_GSM970257_F, xlim=c(-3000, 3000), color="red")
dev.copy(png,'heatmap_TAL1_F.png')
dev.off()

## get binding profile
plotAvgProf(tagMatrix_TAL1_A, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'avgprof_TAL18_A.png')
dev.off()
plotAvgProf(tagMatrix_TAL1_F, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'avgprof_TAL1_F.png')
dev.off()

#peak annotation
plotAnnoPie(peakAnno_TAL1_A)
dev.copy(png,'annopie_TAL1_A.png')
dev.off()
plotAnnoPie(peakAnno_TAL1_F)
dev.copy(png,'annopie_TAL1_F.png')
dev.off()

#TF binding
plotAnnoBar(peakAnno_TAL1_A)
dev.copy(png,'annobar_TAL1_A.png')
dev.off()
plotAnnoBar(peakAnno_TAL1_F)
dev.copy(png,'annobar_TAL1_F.png')
dev.off()


#distribution of TF binding relative to TSS
plotDistToTSS(peakAnno_TAL1_A)
dev.copy(png,'DistToTSS_TAL1_A.png')
dev.off()
plotDistToTSS(peakAnno_TAL1_F)
dev.copy(png,'DistToTSS_TAL1_F.png')
dev.off()
