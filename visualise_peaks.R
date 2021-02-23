library(trackViewer)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
library(org.Hs.eg.db)
Adult <- importScore('GSM970258_GATA1-A.wig.gz', format='WIG')
Fetal <- importScore('GSM970257_GATA1-F.wig.gz', format='WIG')

gene_list = read.delim("hglft_genome_3d617_4fe3c0.bed")

for (row in 1:nrow(gene_list)) { 
  gene_name <- gene_list[row, "name"]
  chr <- gene_list[row, "chr"]
  start_pos <- gene_list[row, "start"] - 1000
  end_pos <- gene_list[row, "end"] + 1000
  
  gr <- GRanges(chr, IRanges(start_pos, end_pos))
  trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg18.knownGene, org.Hs.eg.db, gr=gr)
  
  entrezID <- get(gene_name, org.Hs.egSYMBOL2EG)
  theTrack <-geneTrack(entrezID, TxDb.Hsapiens.UCSC.hg18.knownGene, type="gene")[[1]]
  
  viewerStyle <- trackViewerStyle()
  setTrackStyleParam(Adult, "color", c("green", "black"))
  setTrackStyleParam(Fetal, "color", c("blue", "black"))
  setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
  vp <- viewTracks(trackList(Adult, Fetal, trs), gr=gr, viewerStyle=viewerStyle, autoOptimizeStyle=TRUE)

  }



