library(chipenrich)

## read in list of ribosomal protein genes
genes = read.delim("ribsome_geneset.txt")

data('locusdef.hg19.nearest_tss', package = 'chipenrich.data')
data('tss.hg19', package = 'chipenrich.data')
adult_peaks = read_bed('GSM970258_GATA1-A_hg19.bed')
fetal_peaks = read_bed('GSM970257_GATA1-F_hg19.bed')
adult_assigned_peaks <- assign_peaks(peaks = adult_peaks, locusdef = locusdef.hg19.nearest_tss, tss = tss.hg19)
fetal_assigned_peaks <- assign_peaks(peaks = fetal_peaks, locusdef = locusdef.hg19.nearest_tss, tss = tss.hg19)


adult_RP_peaks <- adult_assigned_peaks[adult_assigned_peaks$nearest_tss_gene_id %in% genes$gene_id, ]
filename <- file.path(getwd(), "RPgenePeaks_GSM970258_A_hg19_chipenrich")
write.csv(x=adult_RP_peaks, file=filename)

fetal_RP_peaks <- fetal_assigned_peaks[fetal_assigned_peaks$nearest_tss_gene_id %in% genes$gene_id, ]
filename <- file.path(getwd(), "RPgenePeaks_GSM970257_F_hg19_chipenrich")
write.csv(x=fetal_RP_peaks, file=filename)

#distribution of peaks
plot_dist_to_tss(peaks = 'GSM970258_GATA1-A_hg19.bed.gz', genome = 'hg19')
dev.copy(png,'adult_dist_to_tss.png')
dev.off()
plot_dist_to_tss(peaks = 'GSM970257_GATA1-F_hg19.bed.gz', genome = 'hg19')
dev.copy(png,'fetal_dist_to_tss.png')
dev.off()

#presence of peaks vs locus length
plot_chipenrich_spline(peaks = 'GSM970258_GATA1-A_hg19.bed.gz', locusdef = 'nearest_tss', genome = 'hg19')
dev.copy(png,'adult_peaks_vs_locuslength.png')
dev.off()
plot_chipenrich_spline(peaks = 'GSM970257_GATA1-F_hg19.bed.gz', locusdef = 'nearest_tss', genome = 'hg19')
dev.copy(png,'fetal_peaks_vs_locuslength.png')
dev.off()

#number of peaks vs locus length
plot_polyenrich_spline(peaks = 'GSM970258_GATA1-A_hg19.bed.gz', locusdef = 'nearest_tss', genome = 'hg19')
dev.copy(png,'adult_num_peaks_vs_locuslength.png')
dev.off()
plot_polyenrich_spline(peaks = 'GSM970257_GATA1-F_hg19.bed.gz', locusdef = 'nearest_tss', genome = 'hg19')
dev.copy(png,'fetal_num_peaks_vs_locuslength.png')
dev.off()

#gene coverage vs locus length
plot_gene_coverage(peaks = 'GSM970258_GATA1-A_hg19.bed.gz', locusdef = 'nearest_tss',  genome = 'hg19')
dev.copy(png,'adult_coverage_vs_locuslength.png')
dev.off()
plot_gene_coverage(peaks = 'GSM970257_GATA1-F_hg19.bed.gz', locusdef = 'nearest_tss',  genome = 'hg19')
dev.copy(png,'fetal_coverage_vs_locuslength.png')
dev.off()