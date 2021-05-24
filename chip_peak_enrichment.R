library(chipenrich)

## read in list of ribosomal protein genes
genes = read.delim("ribsome_geneset.txt")

data('locusdef.hg19.nearest_tss', package = 'chipenrich.data')
data('tss.hg19', package = 'chipenrich.data')
adult_peaks = read_bed('./GSE36994_hg19/GSM970258_GATA1-A_hg19.bed')
fetal_peaks = read_bed('./GSE36994_hg19/GSM970257_GATA1-F_hg19.bed')
adult_assigned_peaks <- assign_peaks(peaks = adult_peaks, locusdef = locusdef.hg19.nearest_tss, tss = tss.hg19)
fetal_assigned_peaks <- assign_peaks(peaks = fetal_peaks, locusdef = locusdef.hg19.nearest_tss, tss = tss.hg19)


adult_RP_peaks <- adult_assigned_peaks[adult_assigned_peaks$nearest_tss_gene_id %in% genes$gene_id, ]
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/RPgenePeaks_GSM970258_A_hg19_chipenrich.csv")
write.csv(x=adult_RP_peaks, file=filename)

fetal_RP_peaks <- fetal_assigned_peaks[fetal_assigned_peaks$nearest_tss_gene_id %in% genes$gene_id, ]
filename <- file.path(getwd(), "/GSE36994_hg19/assigned_peaks/RPgenePeaks_GSM970257_F_hg19_chipenrich.csv")
write.csv(x=fetal_RP_peaks, file=filename)

#distribution of peaks
plot_dist_to_tss(peaks = './GSE36994_hg19/GSM970258_GATA1-A_hg19.bed', genome = 'hg19')
dev.copy(png,'./GSE36994_hg19/graphs/chipenrich_adult_dist_to_tss.png')
dev.off()
plot_dist_to_tss(peaks = './GSE36994_hg19/GSM970257_GATA1-F_hg19.bed', genome = 'hg19')
dev.copy(png,'./GSE36994_hg19/graphs/chipenrich_fetal_dist_to_tss.png')
dev.off()

#presence of peaks vs locus length
plot_chipenrich_spline(peaks = './GSE36994_hg19/GSM970258_GATA1-A_hg19.bed', locusdef = 'nearest_tss', genome = 'hg19')
dev.copy(png,'./GSE36994_hg19/graphs/chipenrich_adult_peaks_vs_locuslength.png')
dev.off()
plot_chipenrich_spline(peaks = './GSE36994_hg19/GSM970257_GATA1-F_hg19.bed', locusdef = 'nearest_tss', genome = 'hg19')
dev.copy(png,'./GSE36994_hg19/graphs/chipenrich_fetal_peaks_vs_locuslength.png')
dev.off()

#number of peaks vs locus length
plot_polyenrich_spline(peaks = './GSE36994_hg19/GSM970258_GATA1-A_hg19.bed', locusdef = 'nearest_tss', genome = 'hg19')
dev.copy(png,'./GSE36994_hg19/graphs/chipenrich_adult_num_peaks_vs_locuslength.png')
dev.off()
plot_polyenrich_spline(peaks = './GSE36994_hg19/GSM970257_GATA1-F_hg19.bed', locusdef = 'nearest_tss', genome = 'hg19')
dev.copy(png,'./GSE36994_hg19/graphs/chipenrich_fetal_num_peaks_vs_locuslength.png')
dev.off()

#gene coverage vs locus length
plot_gene_coverage(peaks = './GSE36994_hg19/GSM970258_GATA1-A_hg19.bed', locusdef = 'nearest_tss',  genome = 'hg19')
dev.copy(png,'./GSE36994_hg19/graphs/chipenrich_adult_coverage_vs_locuslength.png')
dev.off()
plot_gene_coverage(peaks = './GSE36994_hg19/GSM970257_GATA1-F_hg19.bed', locusdef = 'nearest_tss',  genome = 'hg19')
dev.copy(png,'./GSE36994_hg19/graphs/chipenrich_fetal_coverage_vs_locuslength.png')
dev.off()



#### histone marks ####
PolII <- read_bed('./GSE36994_hg19/GSM908069_PolII-A_peaks_hg19.bed')
H3K27ac <- read_bed('./GSE36994_hg19/GSM908051_H3K27ac-A_peaks_hg19.bed')
H3K36me3 <- read_bed('./GSE36994_hg19/GSM908047_H3K36me3-A_peaks_hg19.bed')
H3K27me3 <- read_bed('./GSE36994_hg19/GSM908043_H3K27me3-A_peaks_hg19.bed')
H3K9me3 <- read_bed('./GSE36994_hg19/GSM908041_H3K9me3-A_peaks_hg19.bed')
H3K4me3 <- read_bed('./GSE36994_hg19/GSM908039_H3K4me3-A_peaks_hg19.bed')
H3K4me1 <- read_bed('./GSE36994_hg19/GSM908035_H3K4me1-A_peaks_hg19.bed')

assigned_peaks <- assign_peaks(peaks = PolII, locusdef = locusdef.hg19.nearest_tss, tss = tss.hg19)
RP_peaks <- assigned_peaks[assigned_peaks$nearest_tss_gene_id %in% genes$gene_id, ]
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/PolII_chipenrich.csv")
write.csv(x=adult_RP_peaks, file=filename)

assigned_peaks <- assign_peaks(peaks = H3K27ac, locusdef = locusdef.hg19.nearest_tss, tss = tss.hg19)
RP_peaks <- assigned_peaks[assigned_peaks$nearest_tss_gene_id %in% genes$gene_id, ]
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/H3K27ac_chipenrich.csv")
write.csv(x=adult_RP_peaks, file=filename)

assigned_peaks <- assign_peaks(peaks = H3K36me3, locusdef = locusdef.hg19.nearest_tss, tss = tss.hg19)
RP_peaks <- assigned_peaks[assigned_peaks$nearest_tss_gene_id %in% genes$gene_id, ]
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/H3K36me3_chipenrich.csv")
write.csv(x=adult_RP_peaks, file=filename)

assigned_peaks <- assign_peaks(peaks = H3K27me3, locusdef = locusdef.hg19.nearest_tss, tss = tss.hg19)
RP_peaks <- assigned_peaks[assigned_peaks$nearest_tss_gene_id %in% genes$gene_id, ]
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/H3K27me3_chipenrich.csv")
write.csv(x=adult_RP_peaks, file=filename)

assigned_peaks <- assign_peaks(peaks = H3K4me3, locusdef = locusdef.hg19.nearest_tss, tss = tss.hg19)
RP_peaks <- assigned_peaks[assigned_peaks$nearest_tss_gene_id %in% genes$gene_id, ]
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/H3K4me3_chipenrich.csv")
write.csv(x=adult_RP_peaks, file=filename)

assigned_peaks <- assign_peaks(peaks = H3K9me3, locusdef = locusdef.hg19.nearest_tss, tss = tss.hg19)
RP_peaks <- assigned_peaks[assigned_peaks$nearest_tss_gene_id %in% genes$gene_id, ]
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/H3K9me3_chipenrich.csv")
write.csv(x=adult_RP_peaks, file=filename)

assigned_peaks <- assign_peaks(peaks = H3K4me1, locusdef = locusdef.hg19.nearest_tss, tss = tss.hg19)
RP_peaks <- assigned_peaks[assigned_peaks$nearest_tss_gene_id %in% genes$gene_id, ]
filename <- file.path(getwd(), "./GSE36994_hg19/assigned_peaks/histones/H3K4me1_chipenrich.csv")
write.csv(x=adult_RP_peaks, file=filename)