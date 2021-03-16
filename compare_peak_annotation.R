library(VennDiagram)
library(RColorBrewer)

gene_list = read.delim("RPgenes")
chipseeker_adult <- read.csv('RPgenePeaks_GSM970258_A_hg19')
chipseeker_fetal <- read.csv('RPgenePeaks_GSM970257_F_hg19')
chipenrich_adult <- read.csv('RPgenePeaks_GSM970258_A_hg19_chipenrich')
chipenrich_fetal <- read.csv('RPgenePeaks_GSM970257_F_hg19_chipenrich')
chipanno_adult = read.csv("RPgenePeaks_GSM970258_A_hg19_chipanno")
chipanno_adult$SYMBOL <- gene_list$Symbol[match(chipanno_adult$feature, gene_list$Gene.ID)] 
chipanno_fetal = read.csv("RPgenePeaks_GSM970257_F_hg19_chipanno")
chipanno_fetal$SYMBOL <- gene_list$Symbol[match(chipanno_fetal$feature, gene_list$Gene.ID)] 

adult_manual <- read.delim('Adult_filtered_peaks.bed')
fetal_manual <- read.delim('Fetal_filtered_peaks.bed')

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(x = list(unique(chipseeker_adult$SYMBOL), unique(chipenrich_adult$gene_symbol), unique(chipanno_adult$SYMBOL)), category.names = c("ChIPseeker", "ChIPenrich", "ChIPpeakAnno"), filename = 'gene_overlap_adult.png',
             output=TRUE, imagetype="png" ,
             height = 600 , 
             width = 600 , 
             resolution = 300,
             compression = "lzw",
             
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = myCol,
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.6,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1)



venn.diagram(x = list(chipseeker_fetal$SYMBOL, chipenrich_fetal$gene_symbol, chipanno_fetal$SYMBOL), category.names = c("ChIPseeker", "ChIPenrich", "ChIPpeakAnno"), filename = 'gene_overlap_fetal.png',
             output=TRUE, imagetype="png" ,
             height = 600 , 
             width = 600 , 
             resolution = 300,
             compression = "lzw",
             
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = myCol,
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.6,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1)

adult_genes <- unique(c(chipseeker_adult$SYMBOL, chipenrich_adult$gene_symbol, chipanno_adult$SYMBOL))
fetal_genes <- unique(c(chipseeker_fetal$SYMBOL, chipenrich_fetal$gene_symbol, chipanno_fetal$SYMBOL))
setdiff(adult_genes, fetal_genes)
setdiff(fetal_genes, adult_genes)

venn.diagram(x = list(adult_genes, fetal_genes), category.names = c("Adult", "Fetal"), filename = 'gene_overlap_adult_vs_fetal.png',
             output=TRUE, imagetype="png")

venn.diagram(x = list(unique(adult_manual$gene), fetal_manual$gene), category.names = c("Adult", "Fetal"), filename = 'gene_overlap_manual_filter.png',
             output=TRUE, imagetype="png")

venn.diagram(x = list(adult_genes, unique(adult_manual$gene)), category.names = c("Tools", "Manual"), filename = 'gene_overlap_adult_manual.png',
             output=TRUE, imagetype="png"
)
venn.diagram(x = list(fetal_genes, unique(fetal_manual$gene)), category.names = c("Tools", "Manual"), filename = 'gene_overlap_fetal_manual.png',
             output=TRUE, imagetype="png"
)


annotated_gene_list <- read.delim("RPgenes_2017report")
venn.diagram(x = list(adult_genes, fetal_genes, unique(annotated_gene_list$gene_name)), category.names = c("Adult", "Fetal", "2017 report"), filename = 'gene_overlap_with_2017report.png',
             output=TRUE, imagetype="png" ,
             height = 600 , 
             width = 600 , 
             resolution = 300,
             compression = "lzw",
             
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = myCol,
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.6,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1)

my_annotated_list <- unique(c(adult_genes, fetal_genes))
reported_list <- unique(annotated_gene_list$gene_name)
setdiff(reported_list,my_annotated_list)
setdiff(my_annotated_list, reported_list)

encode_adult <- read.csv('RPgenePeaks.csv')
encode_fetal <- read.csv('RPgenePeaksFetal.csv')
venn.diagram(x = list(unique(encode_adult$SYMBOL), unique(encode_fetal$SYMBOL), unique(annotated_gene_list$gene_name)), category.names = c("ENCODE Adult", "ENCODE Fetal", "2017 report"), filename = 'gene_overlap_ENCODE_2017report.png',
             output=TRUE, imagetype="png" ,
             height = 600 , 
             width = 600 , 
             resolution = 300,
             compression = "lzw",
             
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = myCol,
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.6,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1)


#### histone marks ####
GATA1 <- unique(c(chipseeker_adult$SYMBOL, chipenrich_adult$gene_symbol, chipanno_adult$SYMBOL))
GATA1_manual <- unique(c(GATA1, adult_manual$gene))

#### Pol II and GATA1 ####
chipseeker_PolII <- read.delim('PolII.txt', sep=' ')
chipanno_PolII <- read.csv('PolII_peakanno')
chipenrich_PolII <- read.csv('PolII_chipenrich')
PolII <- unique(c(chipseeker_PolII$SYMBOL, chipanno_PolII$SYMBOL, chipenrich_PolII$gene_symbol))

venn.diagram(x = list(GATA1, PolII), category.names = c("GATA1", "Pol II"), filename = 'gene_overlap_GATA1_PolII.png',
             output=TRUE, imagetype="png")


venn.diagram(x = list(GATA1_manual, PolII), category.names = c("GATA1", "Pol II"), filename = 'gene_overlap_manualGATA1_PolII.png',
             output=TRUE, imagetype="png")

#### enhancer marks H3K4me1 and H3K27ac ####
chipseeker_H3K4me1 <- read.delim('H3K4me1.txt', sep=' ')
chipanno_H3K4me1 <- read.csv('H3K4me1_peakanno')
chipenrich_H3K4me1 <- read.csv('H3K4me1_chipenrich')
H3K4me1 <- unique(c(chipseeker_H3K4me1$SYMBOL, chipanno_H3K4me1$SYMBOL, chipenrich_H3K4me1$gene_symbol))

chipseeker_H3K27ac <- read.delim('H3K27ac.txt', sep=' ')
chipanno_H3K27ac <- read.csv('H3K27ac_peakanno')
chipenrich_H3K27ac <- read.csv('H3K27ac_chipenrich')
H3K27ac <- unique(c(chipseeker_H3K27ac$SYMBOL, chipanno_H3K27ac$SYMBOL, chipenrich_H3K27ac$gene_symbol))

venn.diagram(x = list(GATA1_manual, H3K4me1, H3K27ac), category.names = c("GATA1", "H3K4me1", "H3K27ac"), filename = 'gene_overlap_enhancermarks.png',
             output=TRUE, imagetype="png")

#### promoter marks H3K4me3 ####
chipseeker_H3K4me3 <- read.delim('H3K4me3.txt', sep=' ')
chipanno_H3K4me3 <- read.csv('H3K4me3_peakanno')
chipenrich_H3K4me3 <- read.csv('H3K4me3_chipenrich')
H3K4me3 <- unique(c(chipseeker_H3K4me3$SYMBOL, chipanno_H3K4me3$SYMBOL, chipenrich_H3K4me3$gene_symbol))

venn.diagram(x = list(GATA1_manual, H3K4me3), category.names = c("GATA1", "H3K4me3"), filename = 'gene_overlap_manualGATA1_H3K4me3.png',
             output=TRUE, imagetype="png")

#### transcribed regions in gene bodies H3K36me3 ####
chipseeker_H3K36me3 <- read.delim('H3K36me3.txt', sep=' ')
chipanno_H3K36me3 <- read.csv('H3K36me3_peakanno')
chipenrich_H3K36me3 <- read.csv('H3K36me3_chipenrich')
H3K36me3 <- unique(c(chipseeker_H3K36me3$SYMBOL, chipanno_H3K36me3$SYMBOL, chipenrich_H3K36me3$gene_symbol))

venn.diagram(x = list(GATA1_manual, H3K36me3), category.names = c("GATA1", "H3K36me3"), filename = 'gene_overlap_manualGATA1_H3K36me3.png',
             output=TRUE, imagetype="png")

#### polycomb repression H3K27me3 ####
chipseeker_H3K27me3 <- read.delim('H3K27me3.txt', sep=' ')
chipanno_H3K27me3 <- read.csv('H3K27me3_peakanno')
chipenrich_H3K27me3 <- read.csv('H3K27me3_chipenrich')
H3K27me3 <- unique(c(chipseeker_H3K27me3$SYMBOL, chipanno_H3K27me3$SYMBOL, chipenrich_H3K27me3$gene_symbol))

venn.diagram(x = list(GATA1_manual, H3K27me3), category.names = c("GATA1", "H3K27me3"), filename = 'gene_overlap_manualGATA1_H3K27me3.png',
             output=TRUE, imagetype="png")

#### heterochromatin H3K9me3 ####
chipseeker_H3K9me3 <- read.delim('H3K9me3.txt', sep=' ')
chipanno_H3K9me3 <- read.csv('H3K9me3_peakanno')
chipenrich_H3K9me3 <- read.csv('H3K9me3_chipenrich')
H3K9me3 <- unique(c(chipseeker_H3K9me3$SYMBOL, chipanno_H3K9me3$SYMBOL, chipenrich_H3K9me3$gene_symbol))

venn.diagram(x = list(GATA1_manual, H3K9me3), category.names = c("GATA1", "H3K9me3"), filename = 'gene_overlap_manualGATA1_H3K9me3.png',
             output=TRUE, imagetype="png")
