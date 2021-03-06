## This R script uses ChIPseeker package to analyze the SREBP1 ChIP-seq data
## from A549 and MCF7 cell lines load packages and get annotations

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

# TODO(dlroxe): This is a placeholder 'wrapper' function to enclose
# the code in this file.  Future work should refactor everything it
# contains into more carefully designed functions.
chip_seeker_do_all <- function() {
# tx19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
tx38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- getPromoters(TxDb = tx38, upstream = 3000, downstream = 3000)

# setwd("~/Documents/Bioinformatics_analysis/ChIP analysis/MCF7-SREBP1/MACS")

# TODO(suwu): Update a comment somewhere to explain how big these files are
# and how they were generated; make a plan to get something checked in that
# is capable of generating them (but take care not to exceed any GitHub
# size limits).
files <- list(A549.SREBP1 = "~/Documents/Bioinformatics_analysis/ChIP analysis/A549-SREBP1/MACS/A549-merged-SREBP1-Input_peaks.narrowPeak",
              MCF7.SREBP1 = "~/Documents/Bioinformatics_analysis/ChIP analysis/MCF7-SREBP1/MACS/MCF7-merged-SREBP1-Input_peaks.narrowPeak",
              K562.SREBP1 = "~/Documents/Bioinformatics_analysis/ChIP analysis/MCF7-SREBP1/MACS/K562-merged-SREBP1-Input_peaks.narrowPeak")
print(files)

#############################################################
#####  Profile of ChIP peaks binding to TSS regions  ########
#############################################################
tagMatrixList <- lapply(files, getTagMatrix, windows = promoter)
plotAvgProf(tagMatrixList, xlim = c(-3000, 3000))
tagHeatmap(tagMatrixList, xlim = c(-3000, 3000), color = NULL)

###########################################
#### ChIP peak annotation comparision #####
###########################################
peakAnnoList <- lapply(files, annotatePeak, TxDb = tx38,
                       tssRegion = c(-3000, 3000), verbose = FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

###########################################
##### Functional profiles comparison ######
###########################################
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compREA <- compareCluster(geneCluster    = genes,
                           fun           = "enrichPathway",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
plot(compREA,  showCategory = 20, title  = "Pathway Enrichment Analysis")
compKEGG <- compareCluster(genes,
                           fun          = "enrichKEGG",
                           organism     = "hsa",
                           pvalueCutoff = 0.05)
plot(compKEGG, showCategory = 20, title = "Pathway Enrichment Analysis")

CompareGO_BP <- compareCluster(genes,
                               fun           = "enrichGO",
                               pvalueCutoff  = 0.05,
                               pAdjustMethod = "BH",
                               OrgDb         = org.Hs.eg.db,
                               ont           = "BP",
                               readable      = T)

dotplot(CompareGO_BP, showCategory=30, title="GO - Biological Process")


vennplot(genes)
shuffle(p, TxDb = tx38)
enrichPeakOverlap(queryPeak     = files[[1]],
                  targetPeak    = unlist(files[2:4]),
                  TxDb          = tx38,
                  pAdjustMethod = "BH",
                  nShuffle      = 1000,
                  chainFile     = NULL,
                  verbose       = FALSE)

enrichPeakOverlap(queryPeak     = files[[1]],
                  targetPeak    = unlist(files[1:4]),
                  TxDb          = tx38,
                  pAdjustMethod = "BH",
                  nShuffle      = 50,
                  chainFile     = NULL,
                  verbose       = FALSE)
enrichAnnoOverlap(files[[1]],
                  unlist(files[1:3]),
                  TxDb                 = NULL,
                  pAdjustMethod        = "BH",
                  chainFile            = NULL,
                  distanceToTSS_cutoff = NULL)
files <- getSampleFiles()

data(gcSample)
res <- compareCluster(gcSample, fun = "enrichPathway")
plot(res)

}
