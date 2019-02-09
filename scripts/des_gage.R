# TODO(dlroxe): replace 'source' line with 'library'
# line once everything is organized into a proper
# package.  Despite its location under 'scripts',
# this file assumes that its working directory
# is the project root directory (the parent directory
# of both 'scripts' and 'R')
source(file.path('R', 'DESeqandGAGEforPathwayAnalysis.R'))



###################################################################
## 7. Pathway analysis on ChIP-seq and RNA-seq overlapping genes ##
###################################################################
# TODO(dlroxe): this is a temporary wrapper for code that doesn't
# (yet) run from scratch using only files checked into GitHub
doAll2a <- function() {
  # setwd("~/Documents/Su Wu/Documents/Research/Naar Lab/ChIP-seq")
  # TODO(suwu): check this file into project-data/ ...
  # TODO(dlroxe): ...or generate equivalent data programatically
  ChIP_Seq_and_RNA_Seq_overlap <- read_excel("ChIP-Seq and RNA-Seq overlap.xlsx",
                                             col_names = FALSE)
  colnames(ChIP_Seq_and_RNA_Seq_overlap) <- "symbol"
  View(ChIP_Seq_and_RNA_Seq_overlap)
  ChIP_Seq_and_RNA_Seq_overlap <- merge(ressiRNA,
                                        ChIP_Seq_and_RNA_Seq_overlap,
                                        by = "symbol")
  CR.foldchanges <- ChIP_Seq_and_RNA_Seq_overlap$log2FoldChange
  names(CR.foldchanges) <- ChIP_Seq_and_RNA_Seq_overlap$entrez
  head(CR.foldchanges)
  keggres.sigmet.idx <- gage(CR.foldchanges,
                             gsets    = kegg.sigmet.idx,
                             same.dir = TRUE)
  lapply(keggres.sigmet.idx, head,10)
  fc.go.bp.p <- gage(CR.foldchanges,
                     gsets    = go.bp.gs,
                     same.dir = TRUE)
  lapply(fc.go.bp.p, head,20)
  write.table(fc.go.bp.p$greater,
              file = "ChIP_Seq_and_RNA_Seq_overlap.go.bp.p.greater.txt",
              sep = "\t")
  write.table(fc.go.bp.p$less,
              file = "ChIP_Seq_and_RNA_Seq_overlap.go.bp.p.less.txt",
              sep = "\t")
}

des_gage_do_all()
# doAll2a()
