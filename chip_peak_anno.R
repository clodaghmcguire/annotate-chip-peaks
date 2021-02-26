library(ChIPpeakAnno)
genes = read.delim("ribsome_geneset.txt")

adult_sample <- toGRanges('GSM970258_GATA1-A_hg19.bed', format="BED", header=FALSE) 
fetal_sample <- toGRanges('GSM970257_GATA1-F_hg19.bed', format="BED", header=FALSE, skip=3)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")
annoData

adult_annotatedPeak <- annotatePeakInBatch(adult_sample, AnnotationData = annoData)
adult_annotatedPeak[1:3]
pie1(table(as.data.frame(adult_annotatedPeak)$insideFeature))
adult_annotatedPeak_filtered <- adult_annotatedPeak[adult_annotatedPeak$feature %in% genes$gene_id, ]
filename <- file.path(getwd(), "RPgenePeaks_GSM970258_A_hg19_chipanno")
write.csv(x=adult_annotatedPeak_filtered, file=filename)

fetal_annotatedPeak <- annotatePeakInBatch(fetal_sample, AnnotationData = annoData)
fetal_annotatedPeak[1:3]
pie1(table(as.data.frame(fetal_annotatedPeak)$insideFeature))
fetal_annotatedPeak_filtered <- fetal_annotatedPeak[fetal_annotatedPeak$feature %in% genes$gene_id, ]
filename <- file.path(getwd(), "RPgenePeaks_GSM970257_F_hg19_chipanno")
write.csv(x=fetal_annotatedPeak_filtered, file=filename)

library(org.Hs.eg.db)
adult_genes <- addGeneIDs(adult_annotatedPeak_filtered$feature, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')
fetal_genes <- addGeneIDs(fetal_annotatedPeak_filtered$feature, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"), feature_id_type = 'entrez_id')