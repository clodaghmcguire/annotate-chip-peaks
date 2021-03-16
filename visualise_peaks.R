## visualise peaks for hg18 using trackviewer

library(trackViewer)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
library(org.Hs.eg.db)

Adult <- importScore('GSM970258_GATA1-A.wig.gz', format='WIG')
Fetal <- importScore('GSM970257_GATA1-F.wig.gz', format='WIG')

adult_RP_genes <- read.csv("RPgenePeaks_GSM970258_A")
fetal_RP_genes <- read.csv("RPgenePeaks_GSM970257_F")
gene_list = read.delim("hglft_genome_hg18.bed")

for (row in 1:nrow(gene_list)) { 
  gene_name <- gene_list[row, "name"]
  chr <- gene_list[row, "chr"]
  start_pos <- gene_list[row, "start"] - 10000
  end_pos <- gene_list[row, "end"] + 10000
  
  if(gene_name %in% adult_RP_genes$SYMBOL||gene_name %in% fetal_RP_genes$SYMBOL){
    filename <- paste("./genes/",gene_name,".png", sep ="")
    png(file=filename)
    
    gr <- GRanges(chr, IRanges(start_pos, end_pos))
    ids <- getGeneIDsFromTxDb(gr, TxDb.Hsapiens.UCSC.hg18.knownGene)
    
    symbols <- mget(ids, org.Hs.egSYMBOL)
    genes <-geneTrack(ids, TxDb.Hsapiens.UCSC.hg18.knownGene, symbols, asList=FALSE)
    
    optSty <- optimizeStyle(trackList(Adult, Fetal, genes), theme="safe")
    trackListW <- optSty$tracks
    viewerStyleW <- optSty$style
    viewTracks(trackListW, gr=gr, viewerStyle=viewerStyleW)
    dev.off()
  }

  }


## visualise peaks for hg19 using karyoploteR

library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

adult_RP_genes <- read.csv("RPgenePeaks_GSM970258_A_hg19")
fetal_RP_genes <- read.csv("RPgenePeaks_GSM970257_F_hg19")
gene_list = read.delim("hglft_genome_hg19.bed")

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 15
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

for (row in 1:nrow(gene_list)) { 
  gene_name <- gene_list[row, "name"]
  chr <- gene_list[row, "chr"]
  start_pos <- gene_list[row, "start"] - 10000
  end_pos <- gene_list[row, "end"] + 10000
  
  if(gene_name %in% adult_RP_genes$SYMBOL||gene_name %in% fetal_RP_genes$SYMBOL){
    filename <- paste("./genes/",gene_name,".png", sep ="")
    png(file=filename)

    gene.region <- GRanges(chr, IRanges(start_pos, end_pos))
    kp <- plotKaryotype(zoom = gene.region)
    
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    
    genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                        karyoplot=kp,
                                        plot.transcripts = TRUE, 
                                        plot.transcripts.structure = TRUE)
    
    kp <- plotKaryotype(zoom = gene.region)
    kpPlotGenes(kp, data=genes.data)
    genes.data <- addGeneNames(genes.data)
    genes.data <- mergeTranscripts(genes.data)
    
    kp <- plotKaryotype(zoom = gene.region, cex=1.2, plot.params = pp)
    kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 1.2)
    kpAddBaseNumbers(kp, tick.dist = 10000, minor.tick.dist = 2000,
                     add.units = TRUE, cex=1.3, digits = 6)
    
    total.tracks <- 2
    out.at <- autotrack(1:2, 2, margin = 0.3, r0=0.23)
    
    Adult.bw <- "GSM970258_GATA1-A_hg19.bw"
    at <- autotrack(1, 2, r0=0.23, r1=1, margin = 0.1)
    kp <- kpPlotBigWig(kp, data=Adult.bw, ymax="visible.region", r0=0.35, r1=0.65, col = "cadetblue2")
    computed.ymax <- kp$latest.plot$computed.values$ymax
    kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.35, r1=0.65)
    kpAddLabels(kp, labels = "Adult", r0=0.1, r1=0.65, cex=1.2, label.margin = 0.015)
    
    Fetal.bw <- "GSM970257_GATA1-F_hg19.bw"
    at <- autotrack(2, 2, r0=0.23, r1=1, margin = 0.1)
    kp <- kpPlotBigWig(kp, data=Fetal.bw, ymax="visible.region", r0=0.7, r1=1, col = "darkolivegreen1")
    computed.ymax <- kp$latest.plot$computed.values$ymax
    kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.7, r1=1)
    kpAddLabels(kp, labels = "Fetal", r0=0.6, r1=1, cex=1.2, label.margin = 0.015)
    dev.off()
  }
  
}



## visualise GATA1 and histone marks for hg19 using karyoploteR

library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

adult_RP_genes <- read.csv("RPgenePeaks_GSM970258_A_hg19") # output chipseeker
filtered_bed <- read.delim("Adult_filtered_peaks.bed") # output from manually filtering bed file to within 10kb of RP genes
adult_chipanno <- read.csv("RPgenePeaks_GSM970258_A_hg19_chipanno") # output from chipanno
gene_list = read.delim("RPgenes")
adult_chipanno$SYMBOL <- gene_list$Symbol[match(adult_chipanno$feature, gene_list$Gene.ID)] #add gene symbol to chipanno output
adult_chipenrich <- read.csv("RPgenePeaks_GSM970258_A_hg19_chipenrich") # output from chipenrich
RPgenes <- unique(c(adult_RP_genes$SYMBOL, filtered_bed$gene, adult_chipenrich$gene_symbol, adult_chipanno$SYMBOL))

tools <- unique(c(adult_RP_genes$SYMBOL, adult_chipenrich$gene_symbol, adult_chipanno$SYMBOL))
manual <- unique(filtered_bed$gene)


RPgene_positions = read.delim("hglft_genome_hg19.bed")

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 15
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10
pp$data1outmargin <- 0

for (row in 1:nrow(RPgene_positions)) { 
  gene_name <- RPgene_positions[row, "name"]
  chr <- RPgene_positions[row, "chr"]
  start_pos <- RPgene_positions[row, "start"] - 5000
  end_pos <- RPgene_positions[row, "end"] + 5000
  
  if(gene_name %in% RPgenes){
    filename <- paste("./histone_marks/",gene_name,".png", sep ="")
    png(file=filename, width = 1000, height = 600,
        units = "px", pointsize = 12, bg = "white", res = NA)
    
    gene.region <- GRanges(chr, IRanges(start_pos, end_pos))
    kp <- plotKaryotype(zoom = gene.region)
    
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    
    genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                        karyoplot=kp,
                                        plot.transcripts = TRUE, 
                                        plot.transcripts.structure = TRUE)
    
    genes.data <- addGeneNames(genes.data)
    genes.data <- mergeTranscripts(genes.data)
    
    kp <- plotKaryotype(zoom = gene.region, cex=1.2, plot.params = pp)
    kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 1.2)
    kpAddBaseNumbers(kp, tick.dist = 10000, minor.tick.dist = 2000,
                     add.units = TRUE, cex=1.3, digits = 6)
    
    total.tracks <- 8
    out.at <- autotrack(1:8, 8, margin = 0.3, r0=0.23)
    
    histone.marks <- c(GATA1="GSM970258_GATA1-A_hg19.bw",
                       PolII="GSM908069_PolII-A_hg19.bw",
                       H3K4me1="GSM908035_H3K4me1-A_hg19.bw",
                       H3K27ac="GSM908051_H3K27ac-A_hg19.bw",
                       H3K4me3="GSM908039_H3K4me3-A_hg19.bw",
                       H3K36me3="GSM908047_H3K36me3-A_hg19.bw",
                       H3K27me3="GSM908043_H3K27me3-A_hg19.bw",
                       H3K9me3="GSM908041_H3K9me3-A_hg19.bw"
                       )
    
    colours <- c("#AAAAFF",
                 "#AAFFAA",
                 "#AAFFAA",
                 "#AAFFAA",
                 "#AAFFAA",
                 "#AAFFAA",
                 "#FFAAAA",
                 "#FFAAAA")
    
    for(i in seq_len(length(histone.marks))) {
      bigwig.file <- histone.marks[i]
      at <- autotrack(i, length(histone.marks), r0=0.35, r1=1, margin = 0.5)
      kp <- kpPlotBigWig(kp, data=bigwig.file, ymax="visible.region", 
                         r0=at$r0, r1=at$r1, col = colours[i])
      computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
      kpAxis(kp, ymin=0, ymax=computed.ymax, numticks = 2, r0=at$r0, r1=at$r1)
      kpAddLabels(kp, labels = names(histone.marks)[i], r0=at$r0, r1=at$r1, 
                  cex=1.6, label.margin = 0.035)
    }
    
    dev.off()
  }
  
}


