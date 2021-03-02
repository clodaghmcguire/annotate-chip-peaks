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

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(x = list(chipseeker_adult$SYMBOL, chipenrich_adult$gene_symbol, chipanno_adult$SYMBOL), category.names = c("ChIPseeker", "ChIPenrich", "ChIPpeakAnno"), filename = 'gene_overlap_adult.png',
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