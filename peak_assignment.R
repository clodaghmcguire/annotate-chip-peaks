library(chipenrich)
library(ChIPseeker)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)
library(tidyr)
library(dplyr)
library(VennDiagram)
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
seqlevelsStyle(edb) <- "UCSC"

## read in list of ribosomal protein genes
genes = read.delim("ribsome_geneset.txt")

## datasets can be downloaded from GSE36994
## original FASTQ aligned to hg19 and peaks called using macs3
## annotate peaks
peakAnno_GATA1 <- annotatePeak('./GSE36994_hg19_macs3/macs3_q1e-2_adult/GATA1_macs3_peaks_enrichment5.bed', 
                               tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Hs.eg.db", 
                               addFlankGeneInfo = TRUE, flankDistance = 5000)

## get list of assigned peaks
all_peaks_GATA1 <- as.GRanges(peakAnno_GATA1) #use GRanges rather than data.frame so that after filtering you can use covplot
head(all_peaks_GATA1)

## filter annoPeaks to just ribosomal protein genes
filtered_peaks_GATA1 <- all_peaks_GATA1[all_peaks_GATA1$geneId %in% genes$gene_id, ]
head(filtered_peaks_GATA1)
filename <- file.path(getwd(), "./GSE36994_hg19_macs3/assigned_peaks/GATA1_chipseeker_flank_info_txdb.csv")
write.csv(x=filtered_peaks_GATA1, file=filename, row.names = FALSE)

## filter annoPeaks to get peaks where RP gene is a flanking gene
flank_df <- all_peaks_GATA1 %>%
  separate_rows(flank_geneIds,)

filtered_flank_df <- flank_df[flank_df$flank_geneIds %in% genes$gene_id, ]

filename <- file.path(getwd(), "./GSE36994_hg19_macs3/macs3_q1e-2_adult/assigned_peaks/GATA1_flank.csv")
write.csv(x=filtered_flank_df, file=filename, row.names = FALSE)


GATA1_peak <- readPeakFile('./GSE36994_hg19_macs3/macs3_q1e-2_adult/GATA1_macs3_peaks_enrichment5.bed')


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

#### adult fetal comparison ####
## annotate peaks
peakAnno_GSM970258_A <- annotatePeak('./GSE36994_hg19_macs3/macs3_q1e-2_adult/GATA1_macs3_peaks_enrichment5.bed', 
                                     tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno_GSM970257_F <- annotatePeak('./GSE36994_hg19_macs3/macs3_q1e-2_fetal/GATA1_macs3_enrichment5.bed', 
                                     tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
## get list of assigned peaks
all_peaks_GSM970258_A <- as.GRanges(peakAnno_GSM970258_A) #use GRanges rather than data.frame so that after filtering you can use covplot
head(all_peaks_GSM970258_A)
all_peaks_GSM970257_F <- as.GRanges(peakAnno_GSM970257_F)
head(all_peaks_GSM970257_F)

## filter annoPeaks to just ribosomal protein genes
filtered_peaks_GSM970258_A <- all_peaks_GSM970258_A[all_peaks_GSM970258_A$geneId %in% genes$gene_id, ]
head(filtered_peaks_GSM970258_A)
filtered_peaks_GSM970257_F <- all_peaks_GSM970257_F[all_peaks_GSM970257_F$geneId %in% genes$gene_id, ]
filename <- file.path(getwd(), "./GSE36994_hg19_macs3/assigned_peaks/RPgenePeaks_GSM970258_A_hg19_chipseeker.csv")
write.csv(x=filtered_peaks_GSM970258_A, file=filename, row.names = FALSE)
filename <- file.path(getwd(), "./GSE36994_hg19_macs3/assigned_peaks/RPgenePeaks_GSM970257_F_hg19.csv")
write.csv(x=filtered_peaks_GSM970257_F, file=filename, row.names = FALSE)

GSM970258_A_peak <- readPeakFile('./GSE36994_hg19_macs3/macs3_q1e-2_adult/GATA1_macs3_peaks_enrichment5.bed')
GSM970257_F_peak <- readPeakFile('./GSE36994_hg19_macs3/macs3_q1e-2_fetal/GATA1_macs3_enrichment5.bed')

covplot(GSM970258_A_peak, weightCol = 'V5')
dev.copy(png,'./GSE36994_hg19_macs3/graphs/peak_coverage_GSM970258_A.png')
dev.off()
covplot(filtered_peaks_GSM970258_A, weightCol = 'V5')
dev.copy(png,'./GSE36994_hg19_macs3/graphs/peak_coverage_GSM970258_A_RPgenes.png')
dev.off()

covplot(GSM970257_F_peak, weightCol = 'V5')
dev.copy(png,'./GSE36994_hg19_macs3/graphs/peak_coverage_GSM970257_F.png')
dev.off()
covplot(filtered_peaks_GSM970257_F, weightCol = 'V5')
dev.copy(png,'./GSE36994_hg19_macs3/graphs/peak_coverage_GSM970257_F_RPgenes.png')
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_GSM970258_A <- getTagMatrix(GSM970258_A_peak, windows=promoter)
tagMatrix_GSM970257_F <- getTagMatrix(GSM970257_F_peak, windows=promoter)

## get binding profile
plotAvgProf(tagMatrix_GSM970258_A, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./GSE36994_hg19_macs3/graphs/avgprof_GSM970258_A.png')
dev.off()
plotAvgProf(tagMatrix_GSM970257_F, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./GSE36994_hg19_macs3/graphs/avgprof_GSM970257_F.png')
dev.off()

#peak annotation
plotAnnoPie(peakAnno_GSM970258_A)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/annopie_GSM970258_A.png')
dev.off()
plotAnnoPie(peakAnno_GSM970257_F)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/annopie_GSM970257_F.png')
dev.off()


#TF binding
plotAnnoBar(peakAnno_GSM970258_A)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/annobar_GSM970258_A.png', width=480, height=240)
dev.off()
plotAnnoBar(peakAnno_GSM970257_F)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/annobar_GSM970257_F.png', width=480, height=240)
dev.off()


#distribution of TF binding relative to TSS
plotDistToTSS(peakAnno_GSM970258_A, title = "Distribution of TF binding relative to TSS in adult cells")
dev.copy(png,'./GSE36994_hg19_macs3/graphs/DistToTSS_GSM970258_A.png', width=480, height=240)
dev.off()
plotDistToTSS(peakAnno_GSM970257_F, title = "Distribution of TF binding relative to TSS in fetal cells")
dev.copy(png,'./GSE36994_hg19_macs3/graphs/DistToTSS_GSM970257_F.png', width=480, height=240)
dev.off()

##comparison of adult vs fetal GATA1
GATA1_files <- list(Adult = './GSE36994_hg19_macs3/macs3_q1e-2_adult/GATA1_macs3_peaks_enrichment5.bed', 
                    Fetal = './GSE36994_hg19_macs3/macs3_q1e-2_fetal/GATA1_macs3_enrichment5.bed')

tagMatrixList_GATA1 <- lapply(GATA1_files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList_GATA1, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./GSE36994_hg19_macs3/graphs/avgprof_GATA1.png')
dev.off()
PlotAvgProf(tagMatrixList_GATA1, xlim=c(-3000, 3000), conf=0.95, resample=500, facet="row")
dev.copy(png,'./GSE36994_hg19_macs3/graphs/avgprof_conf_GATA1.png')
dev.off()
peakAnnoList_GATA1 <- lapply(GATA1_files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList_GATA1)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/annobar_GATA1.png', width=480, height=240)
dev.off()
plotDistToTSS(peakAnnoList_GATA1)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/DistToTSS_GATA1.png', width=480, height=240)
dev.off()

GATA1_all_genes <- lapply(peakAnnoList_GATA1, function(i) unique(as.data.frame(i)$geneId))
vennplot(GATA1_all_genes)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/overlap_all_genes.png', width=480, height=480)
dev.off()

GATA1_rp_genes <- lapply(peakAnnoList_GATA1, function(i) unique(as.data.frame(i)$geneId[as.data.frame(i)$geneId %in% genes$gene_id]))
vennplot(GATA1_rp_genes)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/overlap_RP_genes.png', width=480, height=480)
dev.off()

#comparison of adult vs fetal for PolII
PolII_files <- list(PolII_adult = './GSE36994_hg19_macs3/macs3_PolII_adult/PolII_macs3_peaks.narrowPeak', 
                    PolII_fetal = './GSE36994_hg19_macs3/macs3_PolII_fetal/PolII_macs3_peaks.narrowPeak')

tagMatrixList_PolII <- lapply(PolII_files, getTagMatrix, windows=promoter)

plotAvgProf(tagMatrixList_PolII, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.copy(png,'./GSE36994_hg19_macs3/graphs/avgprof_PolII.png')
dev.off()
PlotAvgProf(tagMatrixList_PolII, xlim=c(-3000, 3000), conf=0.95, resample=500, facet="row")
dev.copy(png,'./GSE36994_hg19_macs3/graphs/avgprof_conf_PolII.png')
dev.off()

peakAnnoList_PolII <- lapply(PolII_files, annotatePeak, TxDb=edb, annoDb="org.Hs.eg.db", 
                             tssRegion=c(-1000, 1000), verbose=FALSE)

plotAnnoBar(peakAnnoList_PolII)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/annobar_PolII.png', width=480, height=240)
dev.off()
plotDistToTSS(peakAnnoList_PolII)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/DistToTSS_PolII.png', width=480, height=240)
dev.off()

PolII_rp_gene <- lapply(peakAnnoList_PolII, function(i) unique(as.data.frame(i)$geneId[as.data.frame(i)$geneId %in% genes$gene_id]))
vennplot(PolII_rp_gene)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/overlap_PolII_genes.png', width=480, height=480)
dev.off()

PolII_df <- lapply(peakAnnoList_PolII, function(i) as.data.frame(i)[as.data.frame(i)$ENTREZID %in% genes$gene_id, ])
sapply(names(PolII_df), function (x) write.csv(PolII_df[[x]], file=paste(paste("./GSE36994_hg19_macs/assigned_peaks/", x), "csv", sep="."), row.names = FALSE )   )

plotDistToTSS(list(GATA1 = peakAnno_GSM970258_A, Pol_II = peakAnnoList_PolII$Adult))
dev.copy(png,'./GSE36994_hg19_macs3/graphs/DistToTSS_PolII_GATA1_adult.png', width=480, height=240)
dev.off()
plotDistToTSS(list(GATA1 = peakAnno_GSM970257_F, Pol_II = peakAnnoList_PolII$Fetal))
dev.copy(png,'./GSE36994_hg19_macs3/graphs/DistToTSS_PolII_GATA1_fetal.png', width=480, height=240)
dev.off()
### tool comparison ###
#chipenrich
data('locusdef.hg19.nearest_tss', package = 'chipenrich.data')
data('tss.hg19', package = 'chipenrich.data')
#dropped all columns except position & peak name for chipenrich or read_bed fails
adult_peaks = read_bed('./GSE36994_hg19_macs3/macs3_q1e-2_adult/GATA1_macs3_peaks_enrichment5_forchipenrich.bed') 
fetal_peaks = read_bed('./GSE36994_hg19_macs3/macs3_q1e-2_fetal/GATA1_macs3_enrichment5_forchipenrich.bed')
adult_assigned_peaks <- assign_peaks(peaks = adult_peaks, locusdef = locusdef.hg19.nearest_tss, tss = tss.hg19)
fetal_assigned_peaks <- assign_peaks(peaks = fetal_peaks, locusdef = locusdef.hg19.nearest_tss, tss = tss.hg19)


adult_RP_peaks <- adult_assigned_peaks[adult_assigned_peaks$nearest_tss_gene_id %in% genes$gene_id, ]
filename <- file.path(getwd(), "./GSE36994_hg19_macs3/assigned_peaks/RPgenePeaks_A_hg19_chipenrich.csv")
write.csv(x=adult_RP_peaks, file=filename, row.names = FALSE)

fetal_RP_peaks <- fetal_assigned_peaks[fetal_assigned_peaks$nearest_tss_gene_id %in% genes$gene_id, ]
filename <- file.path(getwd(), "/GSE36994_hg19_macs3/assigned_peaks/RPgenePeaks_F_hg19_chipenrich.csv")
write.csv(x=fetal_RP_peaks, file=filename, row.names = FALSE)

#chippeakanno
adult_sample <- toGRanges('./GSE36994_hg19_macs3/macs3_q1e-2_adult/GATA1_macs3_peaks_enrichment5.bed', format="BED", header=FALSE) 
fetal_sample <- toGRanges('./GSE36994_hg19_macs3/macs3_q1e-2_fetal/GATA1_macs3_enrichment5.bed', format="BED", header=FALSE)
annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")

adult_annotatedPeak <- annotatePeakInBatch(adult_sample, AnnotationData = annoData)
adult_annotatedPeak[1:3]
adult_annotatedPeak_filtered <- adult_annotatedPeak[adult_annotatedPeak$feature %in% genes$gene_id, ]
gene_list <- addGeneIDs(adult_annotatedPeak_filtered$feature, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
adult_annotatedPeak_filtered$SYMBOL <- gene_list$symbol[match(adult_annotatedPeak_filtered$feature, gene_list$entrez_id)] 
filename <- file.path(getwd(), "./GSE36994_hg19_macs3/assigned_peaks/RPgenePeaks_A_hg19_chipanno.csv")
write.csv(x=adult_annotatedPeak_filtered, file=filename, row.names = FALSE)

adult_aCR<-assignChromosomeRegion(adult_annotatedPeak, nucleotideLevel=FALSE, 
                                  precedence=c("Promoters", "immediateDownstream", 
                                               "fiveUTRs", "threeUTRs", 
                                               "Exons", "Introns"), 
                                  TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
barplot(adult_aCR$percentage)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/adult_chipanno_annobar_GATA1.png')
dev.off()
pie1(adult_aCR$percentage)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/adult_chipanno_annopie_GATA1.png')
dev.off()

fetal_annotatedPeak <- annotatePeakInBatch(fetal_sample, AnnotationData = annoData)
fetal_annotatedPeak[1:3]
fetal_annotatedPeak_filtered <- fetal_annotatedPeak[fetal_annotatedPeak$feature %in% genes$gene_id, ]
gene_list <- addGeneIDs(fetal_annotatedPeak_filtered$feature, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
fetal_annotatedPeak_filtered$SYMBOL <- gene_list$symbol[match(fetal_annotatedPeak_filtered$feature, gene_list$entrez_id)] 
filename <- file.path(getwd(), "./GSE36994_hg19_macs3/assigned_peaks/RPgenePeaks_F_hg19_chipanno.csv")
write.csv(x=fetal_annotatedPeak_filtered, file=filename, row.names = FALSE)

fetal_aCR<-assignChromosomeRegion(fetal_annotatedPeak, nucleotideLevel=FALSE, 
                              precedence=c("Promoters", "immediateDownstream", 
                                           "fiveUTRs", "threeUTRs", 
                                           "Exons", "Introns"), 
                              TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
barplot(fetal_aCR$percentage)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/fetal_chipanno_annobar_GATA1.png')
dev.off()
pie1(fetal_aCR$percentage)
dev.copy(png,'./GSE36994_hg19_macs3/graphs/fetal_chipanno_annopie_GATA1.png')
dev.off()

venn.diagram(x = list(unique(filtered_peaks_GSM970258_A$geneId), unique(adult_RP_peaks$gene_id), 
                      unique(adult_annotatedPeak_filtered$feature)), category.names = c("Chipseeker", "Chipenrich", "ChipPeakAnno"), 
             filename = './GSE36994_hg19_macs3/tool_comparison/adult_tool_comparison.png',
             output=TRUE, imagetype="png")

venn.diagram(x = list(unique(filtered_peaks_GSM970257_F$geneId), unique(fetal_RP_peaks$gene_id), 
                      unique(fetal_annotatedPeak_filtered$feature)), category.names = 
               c("Chipseeker", "Chipenrich", "ChipPeakAnno"), 
             filename = './GSE36994_hg19_macs3/tool_comparison/fetal_tool_comparison.png',
             output=TRUE, imagetype="png")


### Dataset comparison ###
#GSE43625
GSE43625_files <- list(GATA1 = './GSE43625_hg19/GSM1067274_Erythroid_GATA1_peaks_hg19.bed.gz', 
              TAL1 = './GSE43625_hg19/GSM1067277_Erythroid_TAL1_peaks_hg19.bed.gz',
              KLF1 = './GSE43625_hg19/GSM1067275_Erythroid_KLF1_peaks_hg19.bed.gz',
              NFE2 = './GSE43625_hg19/GSM1067276_Erythroid_NFE2_peaks_hg19.bed.gz'
)

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

tagMatrixList_GSE43625 <- lapply(GSE43625_files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList_GSE43625, xlim=c(-3000, 3000))
dev.copy(png,'./GSE43625_hg19/graphs/avgprof.png', width=600, height=300)
dev.off()
plotAvgProf(tagMatrixList_GSE43625, xlim=c(-3000, 3000), facet="row")
dev.copy(png,'./GSE43625_hg19/graphs/avgprof_row.png', width=600, height=300)
dev.off()

peakAnnoList_GSE43625 <- lapply(GSE43625_files, annotatePeak, TxDb=txdb, annoDb="org.Hs.eg.db",
                             tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList_GSE43625)
dev.copy(png,'./GSE43625_hg19/graphs/annobar.png', width=600, height=300)
dev.off()
plotDistToTSS(peakAnnoList_GSE43625)
dev.copy(png,'./GSE43625_hg19/graphs/DistToTss.png', width=600, height=300)
dev.off()

GSE43625_rp_genes <- lapply(peakAnnoList_GSE43625, function(i) as.data.frame(i)[as.data.frame(i)$geneId %in% genes$gene_id, ])
sapply(names(GSE43625_rp_genes), function (x) write.csv(GSE43625_rp_genes[[x]], file=paste(paste("./GSE43625_hg19/assigned_peaks/", x), "csv", sep="."), row.names = FALSE )   )

GSE43625_rp_genes_unique <- lapply(GSE43625_rp_genes, function(i) unique(as.data.frame(i)$geneId))
vennplot(GSE43625_rp_genes_unique)
dev.copy(png,'./GSE43625_hg19/graphs/gene_overlap.png')
dev.off()

vennplot(list(GATA1 = GSE43625_rp_genes_unique$GATA1, KLF1 = GSE43625_rp_genes_unique$KLF1))
dev.copy(png,'./GSE43625_hg19/graphs/GATA1_KLF1_overlap.png')
dev.off()

#ENCODE
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

#### comparison of all datasets: ENCODE, GSE36994 and GSE43625 ####

venn.diagram(x = list(unique(filtered_peaks_GSM970258_A$SYMBOL), unique(filtered_peaks_ENCODE_adult$SYMBOL), 
                      unique(rp_gene_peaks$GATA1$SYMBOL)), category.names = 
               c("GSE36994", "ENCODE", "GSE43625"), filename = './GSE36994_hg19_macs3/dataset_comparison/three_datasets_overlap.png',
             output=TRUE, imagetype="png"
)

unique(filtered_peaks_ENCODE_adult$SYMBOL)
unique(filtered_peaks_GSM970258_A$SYMBOL)
unique(rp_gene_peaks$GATA1$SYMBOL)

# compare Chipseeker to selecting all peaks within 10kb window
window <- read.delim('./GSE36994_hg19_macs3/adult_windowed.bed', header = FALSE)
venn.diagram(x = list(unique(filtered_peaks_GSM970258_A$SYMBOL), unique(window$V4)), 
             category.names = c("Chipseeker", "10kb window"), 
             filename = './GSE36994_hg19_macs3/tool_comparison/adult_10kbWindow_comparison.png',
             output=TRUE, imagetype="png", height = 3000, width = 3000, ext.text = FALSE)

window_fetal <- read.delim('./GSE36994_hg19_macs3/fetal_windowed.bed', header = FALSE)
venn.diagram(x = list(unique(filtered_peaks_GSM970257_F$SYMBOL), unique(window_fetal$V4)), 
             category.names = c("Chipseeker", "10kb window"), 
             filename = './GSE36994_hg19_macs3/tool_comparison/fetal_10kbWindow_comparison.png',
             output=TRUE, imagetype="png", height = 3000, width = 3000, ext.text = FALSE)


#### histone modifications ####
adult_histones <- list( 
                    TAL1_adult = './GSE36994_hg19/GSM908055_TAL1-A_hg19.bed',
                    H3K27ac_adult = './GSE36994_hg19/GSM908051_H3K27ac-A_peaks_hg19.bed',
                    H3K36me3_adult = './GSE36994_hg19/GSM908047_H3K36me3-A_peaks_hg19.bed',
                    H3K27me3_adult = './GSE36994_hg19/GSM908043_H3K27me3-A_peaks_hg19.bed',
                    H3K9me3_adult = './GSE36994_hg19/GSM908041_H3K9me3-A_peaks_hg19.bed',
                    H3K4me3_adult = './GSE36994_hg19/GSM908039_H3K4me3-A_peaks_hg19.bed',
                    H3K4me1_adult = './GSE36994_hg19/GSM908035_H3K4me1-A_peaks_hg19.bed')

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

tagMatrixList_adult_histones <- lapply(adult_histones, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList_adult_histones, xlim=c(-3000, 3000))
dev.copy(png,'./GSE36994_hg19/graphs/AvgProf_histones_adult.png')
dev.off()

plotAvgProf(tagMatrixList_adult_histones, xlim=c(-3000, 3000), facet="row")
dev.copy(png,'./GSE36994_hg19/graphs/AvgProf_row_histones_adult.png')
dev.off()

par(mar=c(1,1,1,1))
#tagHeatmap(tagMatrixList_adult_histones, xlim=c(-3000, 3000), color=NULL)
peakAnnoList_adult_histones <- lapply(adult_histones, annotatePeak, TxDb=txdb, annoDb="org.Hs.eg.db",
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList_adult_histones)
dev.copy(png,'./GSE36994_hg19/graphs/AnnoBar_histones_adult.png')
dev.off()

plotDistToTSS(list(H3K4me3 = peakAnnoList_adult_histones$H3K4me3_adult, H3K27ac = peakAnnoList_adult_histones$H3K27ac_adult, 
                   H3K4me1 = peakAnnoList_adult_histones$H3K4me1_adult, H3K36me3 = peakAnnoList_adult_histones$H3K36me3_adult,
                   H3K27me3 = peakAnnoList_adult_histones$H3K27me3_adult, H3K9me3 = peakAnnoList_adult_histones$H3K9me3_adult))
dev.copy(png,'./GSE36994_hg19/graphs/DistToTSS_histones_adult.png')
dev.off()

plotDistToTSS(peakAnnoList_adult_histones$TAL1_adult)
dev.copy(png,'./GSE36994_hg19/graphs/DistToTSS_TAL1_adult.png', width=480, height=240)
dev.off()

#all_genes <- lapply(peakAnnoList_adult_histones, function(i) as.data.frame(i)$geneId)
adult_rp_gene_peaks <- lapply(peakAnnoList_adult_histones, function(i) as.data.frame(i)[as.data.frame(i)$geneId %in% genes$gene_id, ])
sapply(names(adult_rp_gene_peaks), function (x) write.csv(adult_rp_gene_peaks[[x]], file=paste(paste("./GSE36994_hg19/assigned_peaks/histones/", x), "csv", sep="."), row.names = FALSE )   )
rp_genes <- lapply(peakAnnoList_adult_histones, function(i) as.data.frame(i)$geneId[as.data.frame(i)$geneId %in% genes$gene_id])

vennplot(list(GATA1 = unique(filtered_peaks_GSM970258_A$SYMBOL), TAL1 = unique(adult_rp_gene_peaks$TAL1_adult$SYMBOL)))
dev.copy(png,'./GSE36994_hg19/graphs/TAL1_GATA1_overlap_adult.png')
dev.off()

fetal_histones <- list( 
  TAL1_fetal = './GSE36994_hg19/GSM908054_TAL1-F_hg19.bed',
  H3K27ac_fetal = './GSE36994_hg19/GSM908050_H3K27ac-F_peaks_hg19.bed',
  H3K36me3_fetal = './GSE36994_hg19/GSM908046_H3K36me3-F_peaks_hg19.bed',
  H3K27me3_fetal = './GSE36994_hg19/GSM908042_H3K27me3-F_peaks_hg19.bed',
  H3K9me3_fetal = './GSE36994_hg19/GSM908040_H3K9me3-F_peaks_hg19.bed',
  H3K4me3_fetal = './GSE36994_hg19/GSM908038_H3K4me3-F_peaks_hg19.bed',
  H3K4me1_fetal = './GSE36994_hg19/GSM908034_H3K4me1-F_peaks_hg19.bed')

tagMatrixList_fetal_histones <- lapply(fetal_histones, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList_fetal_histones, xlim=c(-3000, 3000))
dev.copy(png,'./GSE36994_hg19/graphs/AvgProf_histones_fetal.png')
dev.off()

plotAvgProf(tagMatrixList_fetal_histones, xlim=c(-3000, 3000), facet="row")
dev.copy(png,'./GSE36994_hg19/graphs/AvgProf_row_histones_fetal.png')
dev.off()
#par(mar=c(1,1,1,1))
#tagHeatmap(tagMatrixList_adult_histones, xlim=c(-3000, 3000), color=NULL)
peakAnnoList_fetal_histones <- lapply(fetal_histones, annotatePeak, TxDb=txdb, annoDb="org.Hs.eg.db",
                                      tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList_fetal_histones)
dev.copy(png,'./GSE36994_hg19/graphs/AnnoBar_histones_fetal.png')
dev.off()

plotDistToTSS(peakAnnoList_fetal_histones)
dev.copy(png,'./GSE36994_hg19/graphs/DistToTSS_histones_fetal.png', width=480, height=480)
dev.off()

plotDistToTSS(peakAnnoList_fetal_histones$TAL1_fetal)
dev.copy(png,'./GSE36994_hg19/graphs/DistToTSS_TAL1_fetal.png', width=480, height=240)
dev.off()

#all_genes <- lapply(peakAnnoList_fetal_histones, function(i) as.data.frame(i)$geneId)
fetal_rp_gene_peaks <- lapply(peakAnnoList_fetal_histones, function(i) as.data.frame(i)[as.data.frame(i)$geneId %in% genes$gene_id, ])
sapply(names(fetal_rp_gene_peaks), function (x) write.csv(fetal_rp_gene_peaks[[x]], file=paste(paste("./GSE36994_hg19/assigned_peaks/histones/", x), "csv", sep="."), row.names = FALSE )   )
rp_genes <- lapply(peakAnnoList_fetal_histones, function(i) as.data.frame(i)$geneId[as.data.frame(i)$geneId %in% genes$gene_id])

vennplot(list(GATA1 = unique(filtered_peaks_GSM970257_F$SYMBOL), TAL1 = unique(fetal_rp_gene_peaks$TAL1_fetal$SYMBOL)))
dev.copy(png,'./GSE36994_hg19/graphs/TAL1_GATA1_overlap_fetal.png')
dev.off()

venn.diagram(x = list(unique(fetal_rp_gene_peaks$H3K4me1_fetal$SYMBOL), 
                      unique(adult_rp_gene_peaks$H3K4me1_adult$SYMBOL)), 
             category.names = c("H3K4me1 fetal", "H3K4me1 adult"), 
             filename = './GSE36994_hg19/graphs/H3K4me1_overlap.png',
             output=TRUE, imagetype="png", height = 3000, width = 3000, ext.text = FALSE)

venn.diagram(x = list(unique(fetal_rp_gene_peaks$H3K27ac_fetal$SYMBOL), 
                      unique(adult_rp_gene_peaks$H3K27ac_adult$SYMBOL)), 
             category.names = c("H3K27ac fetal", "H3K27ac adult"), 
             filename = './GSE36994_hg19/graphs/H3K27ac_overlap.png',
             output=TRUE, imagetype="png", height = 3000, width = 3000, ext.text = FALSE)

venn.diagram(x = list(unique(fetal_rp_gene_peaks$H3K4me3_fetal$SYMBOL), 
                      unique(adult_rp_gene_peaks$H3K4me3_adult$SYMBOL)), 
             category.names = c("H3K4me3 fetal", "H3K4me3 adult"), 
             filename = './GSE36994_hg19/graphs/H3K4me3_overlap.png',
             output=TRUE, imagetype="png", height = 3000, width = 3000, ext.text = FALSE)

venn.diagram(x = list(unique(fetal_rp_gene_peaks$H3K36me3_fetal$SYMBOL), 
                      unique(adult_rp_gene_peaks$H3K36me3_adult$SYMBOL)), 
             category.names = c("H3K36me3 fetal", "H3K36me3 adult"), 
             filename = './GSE36994_hg19/graphs/H3K36me3_overlap.png',
             output=TRUE, imagetype="png", height = 3000, width = 3000, ext.text = FALSE)

venn.diagram(x = list(unique(fetal_rp_gene_peaks$H3K9me3_fetal$SYMBOL), 
                      unique(adult_rp_gene_peaks$H3K9me3_adult$SYMBOL)), 
             category.names = c("H3K9me3 fetal", "H3K9me3 adult"), 
             filename = './GSE36994_hg19/graphs/H3K9me3_overlap.png',
             output=TRUE, imagetype="png", height = 3000, width = 3000, ext.text = FALSE)

venn.diagram(x = list(unique(fetal_rp_gene_peaks$H3K27me3_fetal$SYMBOL), 
                      unique(adult_rp_gene_peaks$H3K27me3_adult$SYMBOL)), 
             category.names = c("H3K27me3 fetal", "H3K27me3 adult"), 
             filename = './GSE36994_hg19/graphs/H3K27me3_overlap.png',
             output=TRUE, imagetype="png", height = 3000, width = 3000, ext.text = FALSE)
