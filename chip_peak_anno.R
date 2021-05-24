library(ChIPpeakAnno)
genes = read.delim("ribsome_geneset.txt")

adult_sample <- toGRanges('./GSE36994_hg19/GSM970258_GATA1-A_hg19.bed', format="BED", header=FALSE) 
fetal_sample <- toGRanges('./GSE36994_hg19/GSM970257_GATA1-F_hg19.bed', format="BED", header=FALSE)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
annoData

adult_annotatedPeak <- annotatePeakInBatch(adult_sample, AnnotationData = annoData)
adult_annotatedPeak[1:3]
pie1(table(as.data.frame(adult_annotatedPeak)$insideFeature))
adult_annotatedPeak_filtered <- adult_annotatedPeak[adult_annotatedPeak$feature %in% genes$gene_id, ]
gene_list <- addGeneIDs(adult_annotatedPeak_filtered$feature, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
adult_annotatedPeak_filtered$SYMBOL <- gene_list$symbol[match(adult_annotatedPeak_filtered$feature, gene_list$entrez_id)] 
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/RPgenePeaks_GSM970258_A_hg19_chipanno.csv")
write.csv(x=adult_annotatedPeak_filtered, file=filename)

fetal_annotatedPeak <- annotatePeakInBatch(fetal_sample, AnnotationData = annoData)
fetal_annotatedPeak[1:3]
pie1(table(as.data.frame(fetal_annotatedPeak)$insideFeature))
fetal_annotatedPeak_filtered <- fetal_annotatedPeak[fetal_annotatedPeak$feature %in% genes$gene_id, ]
gene_list <- addGeneIDs(fetal_annotatedPeak_filtered$feature, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
fetal_annotatedPeak_filtered$SYMBOL <- gene_list$symbol[match(fetal_annotatedPeak_filtered$feature, gene_list$entrez_id)] 
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/RPgenePeaks_GSM970257_F_hg19_chipanno.csv")
write.csv(x=fetal_annotatedPeak_filtered, file=filename)


#### histone marks ####

PolII <- toGRanges('./GSE36994_hg19/GSM908069_PolII-A_peaks_hg19.bed', format="BED", header=FALSE)
H3K27ac <- toGRanges('./GSE36994_hg19/GSM908051_H3K27ac-A_peaks_hg19.bed', format="BED", header=FALSE)
H3K36me3 <- toGRanges('./GSE36994_hg19/GSM908047_H3K36me3-A_peaks_hg19.bed', format="BED", header=FALSE)
H3K27me3 <- toGRanges('./GSE36994_hg19/GSM908043_H3K27me3-A_peaks_hg19.bed', format="BED", header=FALSE)
H3K9me3 <- toGRanges('./GSE36994_hg19/GSM908041_H3K9me3-A_peaks_hg19.bed', format="BED", header=FALSE)
H3K4me3 <- toGRanges('./GSE36994_hg19/GSM908039_H3K4me3-A_peaks_hg19.bed', format="BED", header=FALSE)
H3K4me1 <- toGRanges('./GSE36994_hg19/GSM908035_H3K4me1-A_peaks_hg19.bed', format="BED", header=FALSE)



annotatedPeak <- annotatePeakInBatch(PolII, AnnotationData = annoData)
pie1(table(as.data.frame(annotatedPeak)$insideFeature))
annotatedPeak_filtered <- annotatedPeak[annotatedPeak$feature %in% genes$gene_id, ]
gene_list <- addGeneIDs(annotatedPeak_filtered$feature, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
annotatedPeak_filtered$SYMBOL <- gene_list$symbol[match(annotatedPeak_filtered$feature, gene_list$entrez_id)] 
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/PolII_peakanno.csv")
write.csv(x=annotatedPeak_filtered, file=filename)

annotatedPeak <- annotatePeakInBatch(H3K27ac, AnnotationData = annoData)
pie1(table(as.data.frame(annotatedPeak)$insideFeature))
annotatedPeak_filtered <- annotatedPeak[annotatedPeak$feature %in% genes$gene_id, ]
gene_list <- addGeneIDs(annotatedPeak_filtered$feature, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
annotatedPeak_filtered$SYMBOL <- gene_list$symbol[match(annotatedPeak_filtered$feature, gene_list$entrez_id)] 
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/H3K27ac_peakanno.csv")
write.csv(x=annotatedPeak_filtered, file=filename)

annotatedPeak <- annotatePeakInBatch(H3K36me3, AnnotationData = annoData)
pie1(table(as.data.frame(annotatedPeak)$insideFeature))
annotatedPeak_filtered <- annotatedPeak[annotatedPeak$feature %in% genes$gene_id, ]
gene_list <- addGeneIDs(annotatedPeak_filtered$feature, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
annotatedPeak_filtered$SYMBOL <- gene_list$symbol[match(annotatedPeak_filtered$feature, gene_list$entrez_id)] 
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/H3K36me3_peakanno.csv")
write.csv(x=annotatedPeak_filtered, file=filename)

annotatedPeak <- annotatePeakInBatch(H3K27me3, AnnotationData = annoData)
pie1(table(as.data.frame(annotatedPeak)$insideFeature))
annotatedPeak_filtered <- annotatedPeak[annotatedPeak$feature %in% genes$gene_id, ]
gene_list <- addGeneIDs(annotatedPeak_filtered$feature, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
annotatedPeak_filtered$SYMBOL <- gene_list$symbol[match(annotatedPeak_filtered$feature, gene_list$entrez_id)] 
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/H3K27me3_peakanno.csv")
write.csv(x=annotatedPeak_filtered, file=filename)

annotatedPeak <- annotatePeakInBatch(H3K9me3, AnnotationData = annoData)
pie1(table(as.data.frame(annotatedPeak)$insideFeature))
annotatedPeak_filtered <- annotatedPeak[annotatedPeak$feature %in% genes$gene_id, ]
gene_list <- addGeneIDs(annotatedPeak_filtered$feature, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
annotatedPeak_filtered$SYMBOL <- gene_list$symbol[match(annotatedPeak_filtered$feature, gene_list$entrez_id)] 
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/H3K9me3_peakanno.csv")
write.csv(x=annotatedPeak_filtered, file=filename)

annotatedPeak <- annotatePeakInBatch(H3K4me3, AnnotationData = annoData)
pie1(table(as.data.frame(annotatedPeak)$insideFeature))
annotatedPeak_filtered <- annotatedPeak[annotatedPeak$feature %in% genes$gene_id, ]
gene_list <- addGeneIDs(annotatedPeak_filtered$feature, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
annotatedPeak_filtered$SYMBOL <- gene_list$symbol[match(annotatedPeak_filtered$feature, gene_list$entrez_id)] 
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/H3K4me3_peakanno.csv")
write.csv(x=annotatedPeak_filtered, file=filename)

annotatedPeak <- annotatePeakInBatch(H3K4me1, AnnotationData = annoData)
pie1(table(as.data.frame(annotatedPeak)$insideFeature))
annotatedPeak_filtered <- annotatedPeak[annotatedPeak$feature %in% genes$gene_id, ]
gene_list <- addGeneIDs(annotatedPeak_filtered$feature, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
annotatedPeak_filtered$SYMBOL <- gene_list$symbol[match(annotatedPeak_filtered$feature, gene_list$entrez_id)] 
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/H3K4me1_peakanno.csv")
write.csv(x=annotatedPeak_filtered, file=filename)


##################################
#### mouse data ####
##################################
library(ChIPpeakAnno)
mouse_genes = read.delim("./mouse/mouse_RPgenes.tsv")
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
annoData <- toGRanges(TxDb.Mmusculus.UCSC.mm9.knownGene, feature="gene")

mouse_replicate1 <- toGRanges('./mouse/mouse_replicate1.bed', format="BED", header=FALSE)
mouse_replicate2 <- toGRanges('./mouse/mouse_replicate2.bed', format="BED", header=FALSE)

annotatedPeak <- annotatePeakInBatch(mouse_replicate1, AnnotationData = annoData)
pie1(table(as.data.frame(annotatedPeak)$insideFeature))
annotatedPeak_filtered <- annotatedPeak[annotatedPeak$feature %in% mouse_genes$Gene.ID, ]
gene_list <- addGeneIDs(annotatedPeak_filtered$feature, orgAnn="org.Mm.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
annotatedPeak_filtered$SYMBOL <- gene_list$symbol[match(annotatedPeak_filtered$feature, gene_list$entrez_id)] 
filename <- file.path(getwd(), "./mouse/assigned_peaks/mouse_rep1_peakanno.csv")
write.csv(x=annotatedPeak_filtered, file=filename)

annotatedPeak <- annotatePeakInBatch(mouse_replicate2, AnnotationData = annoData)
pie1(table(as.data.frame(annotatedPeak)$insideFeature))
annotatedPeak_filtered <- annotatedPeak[annotatedPeak$feature %in% mouse_genes$Gene.ID, ]
gene_list <- addGeneIDs(annotatedPeak_filtered$feature, orgAnn="org.Mm.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
annotatedPeak_filtered$SYMBOL <- gene_list$symbol[match(annotatedPeak_filtered$feature, gene_list$entrez_id)] 
filename <- file.path(getwd(), "./mouse/assigned_peaks/mouse_rep2_peakanno.csv")
write.csv(x=annotatedPeak_filtered, file=filename)