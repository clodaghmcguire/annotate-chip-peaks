library(trackViewer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(karyoploteR)

gene_list = read.delim("hglft_genome_hg19.bed")

Adult <- read.delim('./GSE36994_hg19_macs3/macs3_q1e-2_adult/GATA1_macs3_peaks_enrichment5.bed', header = FALSE)
Fetal <- read.delim('./GSE36994_hg19_macs3/macs3_q1e-2_fetal/GATA1_macs3_enrichment5.bed', header = FALSE)

### visualise peaks for hg19 using karyoploteR ###
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 15
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

for (row in 1:nrow(gene_list)) {
  gene_name <- gene_list[row, "name"]
  chr <- gene_list[row, "X.chr"]
  start_pos <- gene_list[row, "start"] - 10000
  end_pos <- gene_list[row, "end"] + 10000

    filename <- paste("./karyoploteR/",gene_name,".png", sep ="")
    png(file=filename, width = 720, height = 480,
        units = "px")

    gene.region <- GRanges(seqnames=chr, ranges=IRanges(start_pos, end_pos))
    kp <- plotKaryotype(zoom = gene.region)

    genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                        karyoplot=kp,
                                        plot.transcripts = TRUE,
                                        plot.transcripts.structure = TRUE)

    kp <- plotKaryotype(zoom = gene.region)
    kpPlotGenes(kp, data=genes.data)
    genes.data <- addGeneNames(genes.data)
    genes.data <- mergeTranscripts(genes.data)

    kp <- plotKaryotype(zoom = gene.region, cex=1.2, plot.params = pp)
    kpPlotGenes(kp, data=genes.data, r0=0, r1=0.2, gene.name.cex = 1.0)
    kpAddBaseNumbers(kp, tick.dist = 10000, minor.tick.dist = 2000,
                     add.units = TRUE, cex=1.3, digits = 6)

    total.tracks <- 2
    out.at <- autotrack(1:2, 2, margin = 0.3, r0=0.23)

    Adult.bw <- "./GSE36994_hg19/GSM970258_GATA1-A_hg19.bw"
    at <- autotrack(1, 2, r0=0.23, r1=1, margin = 0.1)
    kp <- kpPlotBigWig(kp, data=Adult.bw, ymax="visible.region", r0=0.3, r1=0.6, col = "cadetblue2")
    computed.ymax <- kp$latest.plot$computed.values$ymax
    kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.3, r1=0.6)
    kpAddLabels(kp, labels = "Adult", r0=0.2, r1=0.3, cex=1.2, label.margin = 0.015)
    kpPlotRegions(kp, data=Adult, r0=0.25, r1=0.27, col="cadetblue2", avoid.overlapping = TRUE)
    
    Fetal.bw <- "./GSE36994_hg19/GSM970257_GATA1-F_hg19.bw"
    at <- autotrack(2, 2, r0=0.23, r1=1, margin = 0.1)
    kp <- kpPlotBigWig(kp, data=Fetal.bw, ymax="visible.region", r0=0.7, r1=1, col = "darkolivegreen1")
    computed.ymax <- kp$latest.plot$computed.values$ymax
    kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.7, r1=1)
    kpAddLabels(kp, labels = "Fetal", r0=0.6, r1=0.7, cex=1.2, label.margin = 0.015)
    kpPlotRegions(kp, data=Fetal, r0=0.65, r1=0.67, col="darkolivegreen1")
    
    dev.off()

}
# 
# 
# 
# ### visualise GATA1 and histone marks for hg19 using karyoploteR ###
# pp <- getDefaultPlotParams(plot.type=1)
# pp$leftmargin <- 0.15
# pp$topmargin <- 15
# pp$bottommargin <- 15
# pp$ideogramheight <- 5
# pp$data1inmargin <- 10
# pp$data1outmargin <- 0
# 
# for (row in 1:nrow(gene_list)) { 
#   gene_name <- gene_list[row, "name"]
#   chr <- gene_list[row, "X.chr"]
#   start_pos <- gene_list[row, "start"] - 5000
#   end_pos <- gene_list[row, "end"] + 5000
#   print(chr)
#   
#   filename <- paste("./histone_marks/",gene_name,".png", sep ="")
#   png(file=filename, width = 1000, height = 600,
#       units = "px")
#   
#   gene.region <- GRanges(chr, IRanges(start_pos, end_pos))
#   kp <- plotKaryotype(zoom = gene.region)
#   
#   genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,
#                                       karyoplot=kp,
#                                       plot.transcripts = TRUE, 
#                                       plot.transcripts.structure = TRUE)
#   
#   genes.data <- addGeneNames(genes.data)
#   genes.data <- mergeTranscripts(genes.data)
#   
#   kp <- plotKaryotype(zoom = gene.region, cex=1.2, plot.params = pp)
#   kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 1.2)
#   kpAddBaseNumbers(kp, tick.dist = 10000, minor.tick.dist = 2000,
#                    add.units = TRUE, cex=1.3, digits = 6)
#   
#   total.tracks <- 9
#   out.at <- autotrack(1:9, 9, margin = 0.3, r0=0.23)
#   
#   histone.marks <- c(GATA1="./GSE36994_hg19/GSM970258_GATA1-A_hg19.bw",
#                      TAL1="./GSE36994_hg19/GSM908055_TAL1-A_hg19.bw",
#                      PolII="./GSE36994_hg19/GSM908069_PolII-A_hg19.bw",
#                      H3K4me1="./GSE36994_hg19/GSM908035_H3K4me1-A_hg19.bw",
#                      H3K27ac="./GSE36994_hg19/GSM908051_H3K27ac-A_hg19.bw",
#                      H3K4me3="./GSE36994_hg19/GSM908039_H3K4me3-A_hg19.bw",
#                      H3K36me3="./GSE36994_hg19/GSM908047_H3K36me3-A_hg19.bw",
#                      H3K27me3="./GSE36994_hg19/GSM908043_H3K27me3-A_hg19.bw",
#                      H3K9me3="./GSE36994_hg19/GSM908041_H3K9me3-A_hg19.bw"
#                      )
#   
#   colours <- c("#AAAAFF",
#                "#AAAAFF",
#                "#AAFFAA",
#                "#AAFFAA",
#                "#AAFFAA",
#                "#AAFFAA",
#                "#AAFFAA",
#                "#FFAAAA",
#                "#FFAAAA")
#   
#   for(i in seq_len(length(histone.marks))) {
#     bigwig.file <- histone.marks[i]
#     at <- autotrack(i, length(histone.marks), r0=0.35, r1=1, margin = 0.5)
#     kp <- kpPlotBigWig(kp, data=bigwig.file, ymax="visible.region", 
#                        r0=at$r0, r1=at$r1, col = colours[i])
#     computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
#     kpAxis(kp, ymin=0, ymax=computed.ymax, numticks = 2, r0=at$r0, r1=at$r1)
#     kpAddLabels(kp, labels = names(histone.marks)[i], r0=at$r0, r1=at$r1, 
#                 cex=1.6, label.margin = 0.035)
#     }
#     
#     dev.off()
#   
# }
# 
# ### visualise peaks for macs3 re-analysis using karyoploteR ###
# pp <- getDefaultPlotParams(plot.type=1)
# pp$leftmargin <- 0.15
# pp$topmargin <- 15
# pp$bottommargin <- 15
# pp$ideogramheight <- 5
# pp$data1inmargin <- 10
# 
# for (row in 1:nrow(gene_list)) { 
#   gene_name <- gene_list[row, "name"]
#   chr <- gene_list[row, "X.chr"]
#   start_pos <- gene_list[row, "start"] - 10000
#   end_pos <- gene_list[row, "end"] + 10000
#   
#   filename <- paste("./GSE36994_hg19_macs3/karyploteR/",gene_name,".png", sep ="")
#   png(file=filename, width = 720, height = 480,
#       units = "px")
#   
#   gene.region <- GRanges(seqnames=chr, ranges=IRanges(start_pos, end_pos))
#   kp <- plotKaryotype(zoom = gene.region)
#   
#   genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,
#                                       karyoplot=kp,
#                                       plot.transcripts = TRUE, 
#                                       plot.transcripts.structure = TRUE)
#   
#   kp <- plotKaryotype(zoom = gene.region)
#   kpPlotGenes(kp, data=genes.data)
#   genes.data <- addGeneNames(genes.data)
#   genes.data <- mergeTranscripts(genes.data)
#   
#   kp <- plotKaryotype(zoom = gene.region, cex=1.2, plot.params = pp)
#   kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 1.2)
#   kpAddBaseNumbers(kp, tick.dist = 10000, minor.tick.dist = 2000,
#                    add.units = TRUE, cex=1.3, digits = 6)
#   
#   total.tracks <- 2
#   out.at <- autotrack(1:2, 2, margin = 0.3, r0=0.23)
#   
#   Adult.bw <- "./GSE36994_hg19_macs3/GSE36994_GATA1.bigwig"
#   at <- autotrack(1, 2, r0=0.23, r1=1, margin = 0.1)
#   kp <- kpPlotBigWig(kp, data=Adult.bw, ymax="visible.region", r0=0.35, r1=0.65, col = "cadetblue2")
#   computed.ymax <- kp$latest.plot$computed.values$ymax
#   kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.35, r1=0.65)
#   kpAddLabels(kp, labels = "GATA1", r0=0.2, r1=0.65, cex=1.2, label.margin = 0.015)
#   dev.off()
#   
# }
