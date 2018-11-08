## This R script uses ChIPseeker package to analyze the SREBP1 ChIP-seq data from A549 and MCF7 cell lines
##load packages and get annotations
source("https://bioconductor.org/biocLite.R")
# biocLite("openssl")
# biocLite("GenomicFeatures")
biocLite("ChIPseeker")
biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
biocLite("clusterProfiler")
biocLite("org.Hs.eg.db")
biocLite("ReactomePA")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

tx19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
tx38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- getPromoters(TxDb=tx38, upstream=3000, downstream=3000)

setwd("~/Documents/Bioinformatics_analysis/ChIP analysis/MCF7-SREBP1/MACS")

files <- list(A549.SREBP1.rep1 = "~/Documents/Bioinformatics_analysis/ChIP analysis/A549-SREBP1/MACS/A549-SREBP1-1-Input-1_peaks.narrowPeak",
              A549.SREBP1.rep2 = "~/Documents/Bioinformatics_analysis/ChIP analysis/A549-SREBP1/MACS/A549-SREBP1-2-Input-1_peaks.narrowPeak",
              MCF7.SREBP1.rep1 = "~/Documents/Bioinformatics_analysis/ChIP analysis/MCF7-SREBP1/MACS/MCF7-SREBP1-1-Input-1_peaks.narrowPeak",
              MCF7.SREBP1.rep2 = "~/Documents/Bioinformatics_analysis/ChIP analysis/MCF7-SREBP1/MACS/MCF7-SREBP1-2-Input-1_peaks.narrowPeak")
 
print(files)


files <- list(A549.SREBP1 = "~/Documents/Bioinformatics analysis/ChIP analysis/A549-SREBP1/MACS/A549-merged-SREBP1-Input_peaks.narrowPeak",
              MCF7.SREBP1 = "~/Documents/Bioinformatics analysis/ChIP analysis/MCF7-SREBP1/MACS/MCF7-merged-SREBP1-Input_peaks.narrowPeak")

#############################################################
#####  Profile of ChIP peaks binding to TSS regions  ########
#############################################################
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

###########################################
#### ChIP peak annotation comparision #####
###########################################
peakAnnoList <- lapply(files, annotatePeak, TxDb=tx38,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

###########################################
##### Functional profiles comparison ######
###########################################
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compREA <- compareCluster(geneCluster   = genes,
                           fun           = "enrichPathway",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
plot(compREA, showCategory = 20, title = "Pathway Enrichment Analysis")
compKEGG <- compareCluster(genes, 
                           fun="enrichKEGG",
                           organism="hsa", 
                           pvalueCutoff=0.05)
plot(compKEGG, showCategory = 20, title = "Pathway Enrichment Analysis")

CompareGO_BP <- compareCluster(genes,
                               fun="enrichGO", 
                               pvalueCutoff=0.05, 
                               pAdjustMethod="BH", 
                               OrgDb=org.Hs.eg.db,
                               ont="BP",
                               readable=T)

dotplot(CompareGO_BP, showCategory=30, title="GO - Biological Process")


vennplot(genes)
shuffle(p, TxDb=tx38)
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
enrichAnnoOverlap(files[[1]], unlist(files[1:3]), TxDb = NULL, pAdjustMethod = "BH",
                  chainFile = NULL, distanceToTSS_cutoff = NULL)
files <- getSampleFiles()

require(clusterProfiler)
data(gcSample)
res <- compareCluster(gcSample, fun="enrichPathway")
plot(res)

############################################################################
setwd("~/Documents/Bioinformatics tools/ChIPseeker/hg19")
hg19 <- getGEOInfo(genome="hg19", simplify=TRUE)
head(hg19)
hg38 <- getGEOInfo(genome="hg38", simplify=TRUE)
head(hg38)
# download BED supplementary files of a list of GSM accession numbers
# GSM733656_hg19_wgEncodeBroadHistoneK562H3k27acStdPk.broadPeak.gz
downloadGSMbedFiles('GSM733656', destDir = "hg19")
K562_hg19_H3k27ac <-"GSM733656_hg19_wgEncodeBroadHistoneK562H3k27acStdPk.broadPeak"
K562_hg19_H3k27ac <- readPeakFile(K562_hg19_H3k27ac)
# GSM733778_hg19_wgEncodeBroadHistoneK562H3k9acStdPk.broadPeak.gz
downloadGSMbedFiles('GSM733778', destDir = "hg19")
K562_hg19_H3k9ac <- "GSM733778_hg19_wgEncodeBroadHistoneK562H3k9acStdPk.broadPeak"
K562_hg19_H3k9ac <- readPeakFile(K562_hg19_H3k9ac)
# GSM1003574_hg19_wgEncodeBroadHistoneK562Cbpsc369Pk.broadPeak.gz
downloadGSMbedFiles('GSM1003574', destDir = "hg19")
K562_hg19_Cbp <- "GSM1003574_hg19_wgEncodeBroadHistoneK562Cbpsc369Pk.broadPeak"
K562_hg19_Cbp <- readPeakFile(K562_hg19_Cbp)
# GSM1003583_hg19_wgEncodeBroadHistoneK562P300StdPk.broadPeak.gz
downloadGSMbedFiles('GSM1003583', destDir = "hg19")
K562_hg19_P300 <- "GSM1003583_hg19_wgEncodeBroadHistoneK562P300StdPk.broadPeak"
K562_hg19_P300 <- readPeakFile(K562_hg19_P300)

# GSM1003578_hg19_wgEncodeBroadHistoneA549H3k27acEtoh02Pk.broadPeak.gz
downloadGSMbedFiles('GSM1003578', destDir = "hg19")
A549_hg19_H3k27ac <- "GSM1003578_hg19_wgEncodeBroadHistoneA549H3k27acEtoh02Pk.broadPeak"
A549_hg19_H3k27ac <- readPeakFile(A549_hg19_H3k27ac)




GSM1003578
Cbp&P300 <- enrichPeakOverlap(hg19_Cbp, hg19_P300, TxDb = NULL, pAdjustMethod = "BH",
                  nShuffle = 1000, chainFile = NULL, pool = TRUE,
                  mc.cores = detectCores() - 1, verbose = TRUE)



######################################################################
K562_hg19_Cbp <- "GSM1003574_hg19_wgEncodeBroadHistoneK562Cbpsc369Pk.broadPeak"
K562_hg19_Cbp <- readPeakFile(K562_hg19_Cbp)
K562_hg19_Cbp.tagMatrix <- getTagMatrix(K562_hg19_Cbp, windows=promoter)
tagHeatmap(K562_hg19_Cbp.tagMatrix, xlim=c(-3000, 3000), color="red")

peak <- K562_hg19_Cbp
peak

peakHeatmap(K562_hg19_Cbp, TxDb=txdb, upstream=3000, downstream=3000, color="red")
plotAvgProf(K562_hg19_Cbp.tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(K562_hg19_Cbp.tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)
peakAnno <- annotatePeak(K562_hg19_Cbp, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")


pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
head(pathway1, 2)
gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)
head(pathway2, 2)
dotplot(pathway2)

