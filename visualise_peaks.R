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



