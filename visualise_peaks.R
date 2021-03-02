## visualise peaks for hg18 using trackviewer

library(trackViewer)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
library(org.Hs.eg.db)

Adult <- importScore('GSM970258_GATA1-A.wig.gz', format='WIG')
Fetal <- importScore('GSM970257_GATA1-F.wig.gz', format='WIG')

adult_RP_genes <- read.csv("RPgenePeaks_GSM970258_A")
fetal_RP_genes <- read.csv("RPgenePeaks_GSM970257_F")
gene_list = read.delim("hglft_genome_3d617_4fe3c0.bed")

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

adult_RP_genes <- read.csv("RPgenePeaks_GSM970258_A")
fetal_RP_genes <- read.csv("RPgenePeaks_GSM970257_F")
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


