# The following script aim to find the differentially expressed
# genes between siSREBF1 and control siRNA treatment in melanoma
# HT-144 cells with DESeq package. And then perform pathway
# analysis on the differentially expressed genes.

library(AnnotationDbi)
library(Biobase)
library(BiocGenerics)
library(BiocStyle)
library(DESeq2)
library(GenomeInfoDb)
library(GenomicRanges)
library(IRanges)
library(Matrix)
library(ReactomePA)
library(RColorBrewer)
library(Rcpp)
library(RcppArmadillo)
library(RODBC)
library(S4Vectors)
library(SummarizedExperiment)
library(ashr)
library(apeglm)
library(biomaRt)
library(colorspace)
library(dplyr)
library(lattice)
library(geneplotter)
library(gage)
library(gageData)
library(genefilter)
library(ggplot2)
library(gplots)
library(org.Hs.eg.db)
library(parallel)
library(pathview)
library(plyr)
library(readxl)
library(rmarkdown)
library(stats4)
library(stringr)
library(survival)

data(kegg.gs)

################################
## 2 Preparing count matrices ##
################################

# TODO(dlroxe): This function is copied from the FactoMineR file.
# Find a single place to define it in common for all callers.
.readGeneCounts <- function () {
  # obtain the count table of the experiment directly from a pre-saved file: gene-counts.csv.
  # The RNA-seq was aligned to human reference genome Hg38 by STAR aligner
  # read processed RNA-seq read data from file testseq.csv.
  testseq <- read.csv(file.path("project-data", "gene-counts-from-Makefile.csv"))
  # Use the column one (Ensemble names) as columnn names.
  testseq <- data.frame(testseq[,-1], row.names=testseq[,1])
  # Remove the first four rows (N_unmapped,N_multimapping,N_noFeature and N_ambiguous)
  testseq <- data.frame(testseq[c(-1,-2,-3,-4),])

  ## remove non-numeric 'symbol col' 25, leaving 4 col X 6 tests
  testseq <- testseq[-25]

  return(testseq)
}


##################################################################################
## 3 Compare siNeg and siBF1 RNA-seq results for differentially expressed genes ##
##################################################################################

# TODO(dlroxe): again, this is exactly the same function as
# found in FactoMineR
.getGuideData <- function(testseq) {
  ## check the distribution of RNA-Seq reads
  par(mar=c(3,12,2,1))
  boxplot(testseq, outline=FALSE, horizontal=TRUE, las=1)
  ## Remove rows in which there are no reads or nearly no reads
  guideData <- testseq[rowSums(testseq)>1,]
  head(guideData)
  dim(guideData)
  ## check how the data distribution with boxplot after removing rows with no read
  par(mar=c(3,12,2,1))
  boxplot(guideData, outline=FALSE, horizontal=TRUE, las=1)
  return(guideData)
}

# TODO(dlroxe): again, same function as other file
.getDDS <- function (guideData, guideDesign, condition) {
  ## Construct DESeqDataSet with the count matrix, countData, and the sample information, colData
  dds <- DESeqDataSetFromMatrix(countData = guideData,
                                colData = guideDesign,
                                design = ~ condition)
  dds
  head(assay(dds))
  return(dds)
}

########################################
## 4 Exploring and exporting results ##
########################################

###############################################################################
## 4.1 standard analysis to make MA-plot from base means and log fold changes##
###############################################################################
.getDDSRes <- function(ddsDEsiRNA) {
  ressiRNA<-results(ddsDEsiRNA) # default alpha = 0.1
  ## p-values and adjusted p-values
  ## We can order our results table by the smallest p value:
  # ressiRNA<-ressiRNA[order(ressiRNA$log2FoldChange),]
  ressiRNA <- ressiRNA[order(ressiRNA$pvalue),]
  ## Information about which variables and tests were used can be found by calling the function mcols on the results object.
  mcols(ressiRNA)$description
  ## We can summarize some basic tallies using the summary function
  summary(ressiRNA)
  ## How many adjusted p-values were less than 0.1?
  sum(ressiRNA$padj < 0.1, na.rm=TRUE)

  ## Note that the results function automatically performs independent filtering based on the mean of normalized counts for each gene,
  ## optimizing the number of genes which will have an adjusted p value below a given FDR cutoff, alpha = 0.05 here
  res05 <- results(ddsDEsiRNA, alpha=0.05)
  summary(res05)
  sum(res05$padj < 0.05, na.rm=TRUE)

  ## plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet.
  ## Points will be colored red if the adjusted p value is less than 0.1.
  ## Points which fall out of the window are plotted as open triangles pointing either up or down.
  plotMA(ressiRNA, main="DESeq2 siRNA", ylim=c(-2,2))
  ## It is more useful visualize the MA-plot for the shrunken log2 fold changes,
  ## which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
  resultsNames(ddsDEsiRNA)
  # normal is the the original DESeq2 shrinkage estimator, an adaptive normal prior
  resLFC <- lfcShrink(ddsDEsiRNA, coef=2)
  # par(mfrow=c(1,3), mar=c(4,4,2,1))
  xlim <- c(1,1e5); ylim <- c(-2.2,2.2)
  plotMA(resLFC,xlab = "mean of normalized counts", ylim=ylim, main="normal",cex=.8)
  # plotMA(resLFC, xlim=xlim, xlab = "mean of normalized counts", ylim=ylim, main="normal")
  return(ressiRNA)
}

# TODO(dlroxe): again, mostly the same as the other file
###########################################################
## 4.2 Add Entrez IDs, gene symbols, and full gene names ##
###########################################################
.addGeneIdentifiers <- function(ressiRNA) {
  ressiRNA <- data.frame(ressiRNA)
  columns(org.Hs.eg.db)
  ressiRNA$symbol = mapIds(org.Hs.eg.db,
                           keys=row.names(ressiRNA),
                           column="SYMBOL",
                           keytype="ENSEMBL",
                           multiVals="first")
  ressiRNA$entrez = mapIds(org.Hs.eg.db,
                           keys=row.names(ressiRNA),
                           column="ENTREZID",
                           keytype="ENSEMBL",
                           multiVals="first")
  ressiRNA$name =   mapIds(org.Hs.eg.db,
                           keys=row.names(ressiRNA),
                           column="GENENAME",
                           keytype="ENSEMBL",
                           multiVals="first")
  summary(ressiRNA)
  head(ressiRNA, 10)
  return(ressiRNA)
}

########################################
## 4.3 Exporting results to CSV files ##
########################################
## Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function,
## followed by the write.csv function.
.exportCSV <- function(ressiRNA) {
  resSig <- subset(ressiRNA, pvalue < 0.05)
  resSig
  write.csv(as.data.frame(resSig),
            file="siBF1_siNeg_sigresults.csv")
  write.csv(as.data.frame(ressiRNA),
            file="siBF1_siNeg_results.csv")
}

##############################################################
## 5. Heatmap for top variable genes across siNeg and siBF1 ##
##############################################################
# The regularized log transform can be obtained using the rlog() function.
# The default “blinds” the normalization to the design.
# topVarGenes looks at the row variance of the transformed values regardless of which samples come from which group.
# Differential expression testing asks whether the difference across group is large relative to the within-group variance.
# So these are different ways of ranking genes.
# calculate the variance for each gene, # select the top 50 genes by variance
.makeHeatMap <- function(ddsDEsiRNA) {
  rldsiRNA <- rlog(ddsDEsiRNA,blind=TRUE)
  vdsiRNA <- varianceStabilizingTransformation(ddsDEsiRNA, blind=TRUE)
  topVarGenes <- head(order(rowVars(assay(vdsiRNA)), decreasing=TRUE), 500)
  mat <- assay(rldsiRNA)[topVarGenes, ]
  mat <- mat - rowMeans(mat)
  heatmap.2(mat,labRow = NA,labCol=NA, scale="row",lhei = c(2, 8),
            trace="none", dendrogram="none", margins=c(2,10),
            col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
            ColSideColors = c(Control="gray", DPN="darkgreen", OHT="orange")[colData(rldsiRNA)$condition ])
}

######################################################
## 6. KEGG pathway analysis on DESeq siNeg vs siBF1 ##
######################################################
.keggAnalysis <- function(ressiRNA) {
  ## Generally Applicable Gene-set/Pathway Analysis (GAGE)
  ## check for coordinated differential expression over gene sets instead of changes of individual genes.
  lapply(kegg.gs[1:3],head)
  ## kegg.sets.hs is a named list of 229 elements
  ## Each element is a character vector of member gene Entrez IDs for a single KEGG pathway
  data(kegg.sets.hs)
  ## sigmet.idx.hs is an index of numbers of sinaling and metabolic pathways in kegg.set.gs.
  data(sigmet.idx.hs)
  ## kegg.sets.hs[sigmet.idx.hs] gives you the “cleaner” gene sets of signalling and metabolic pathways only.
  kegg.sets.hs.sigmet <-  kegg.sets.hs[sigmet.idx.hs]
  head(kegg.sets.hs.sigmet, 3)
  
  ## Generate a named vector of fold changes, where the names of the values are the Entrez gene IDs, for the gage() function
  foldchanges = ressiRNA$log2FoldChange
  names(foldchanges) = ressiRNA$entrez
  head(foldchanges)
  
  # Look at top upregulated (greater) and downregulated (less) genes with statatistics.
  # Get KEGG pathway with both metabolism and signaling pathways
  kg.hsa=kegg.gsets("hsa")
  kegg.sigmet.idx=kg.hsa$kg.sets[kg.hsa$sigmet.idx]
  keggres.sigmet.idx = gage(foldchanges, gsets=kegg.sigmet.idx, same.dir=TRUE)
  lapply(keggres.sigmet.idx, head,20)
  write.table(keggres.sigmet.idx$greater, file = "keggres.sigmet.idx.greater.txt",sep = "\t")
  write.table(keggres.sigmet.idx$less, file = "keggres.sigmet.idx.less.txt",sep = "\t")

  # Get KEGG pathway with only metabolism pathways
  kegg.met.idx=kg.hsa$kg.sets[kg.hsa$met.idx]
  keggres.met.idx = gage(foldchanges, gsets=kegg.met.idx, same.dir=TRUE)
  lapply(keggres.met.idx, head,10)
  # Get KEGG pathway with only signaling pathways
  kegg.sig.idx=kg.hsa$kg.sets[kg.hsa$sig.idx]
  keggres.sig.idx = gage(foldchanges, gsets=kegg.sig.idx, same.dir=TRUE)
  lapply(keggres.sig.idx, head,10)
  return(foldchanges)
}

####################################################
## 6. GO pathway analysis on DESeq siNeg vs siBF1 ##
####################################################
.pathwayAnalysis <- function(foldchanges) {
  # set up GO database
  data(go.gs)
  lapply(go.gs[1:3],head)
  go.hs <- go.gsets(species="human")
  go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP] # “Biological Process”
  go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF] # “Molecular Function”
  go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC] # “Cellular Component”
  # GO pathway analysis with Biological Process
  fc.go.bp.p <- gage(foldchanges, gsets = go.bp.gs,same.dir=TRUE)
  lapply(fc.go.bp.p, head,20)
  write.table(fc.go.bp.p$greater, file = "fc.go.bp.p.greater.txt",sep = "\t")
  write.table(fc.go.bp.p$less, file = "fc.go.bp.p.less.txt",sep = "\t")

  # GO pathway analysis with Molecular Process
  fc.go.mf.p <- gage(foldchanges, gsets = go.mf.gs)
  lapply(fc.go.mf.p, head,10)
  # GO pathway analysis with Cellular Process
  fc.go.cc.p <- gage(foldchanges, gsets = go.cc.gs)
  lapply(fc.go.cc.p, head,10)
}

doAll2 <- function() {
testseq <- .readGeneCounts()
guideDatasiRNA <- .getGuideData(
  data.frame(testseq[,c(5,6,7,8,9,10,11,12)]))

condition <- c(rep("siNeg",4),rep("siBF1",4))
guideDesignsiRNA <- data.frame(row.names = colnames(guideDatasiRNA),
                               condition = condition)

## specifying the reference level:
ddssiRNA <- .getDDS(guideDatasiRNA, guideDesignsiRNA, condition)
ddssiRNA$condition <- relevel(ddssiRNA$condition, ref = "siNeg")
ddssiRNA

ddsDEsiRNA <- DESeq(ddssiRNA)
ressiRNA <- .getDDSRes(ddsDEsiRNA)

ressiRNA <- .addGeneIdentifiers(ressiRNA)

.exportCSV(ressiRNA)

.makeHeatMap(ddsDEsiRNA)

foldchanges <- .keggAnalysis(ressiRNA)

.pathwayAnalysis(foldchanges)

}

# TODO(dlroxe): doAll2a() is a temporary wrapper for code that doesn't
# (yet) run from scratch using only files checked into GitHub
doAll2a <- function() {
########################################################################
###### Pathway analysis on ChIP-seq and RNA-seq overlapping genes ######
########################################################################
#setwd("~/Documents/Su Wu/Documents/Research/Naar Lab/ChIP-seq")
# TODO(suwu): check this file into project-data/ ...
# TODO(dlroxe): ...or generate equivalent data programatically
ChIP_Seq_and_RNA_Seq_overlap <- read_excel("ChIP-Seq and RNA-Seq overlap.xlsx",
                                            col_names = FALSE)
colnames(ChIP_Seq_and_RNA_Seq_overlap) <- "symbol"
View(ChIP_Seq_and_RNA_Seq_overlap)
ChIP_Seq_and_RNA_Seq_overlap <- merge(ressiRNA, ChIP_Seq_and_RNA_Seq_overlap, by="symbol")
CR.foldchanges = -ChIP_Seq_and_RNA_Seq_overlap$log2FoldChange
names(CR.foldchanges) = ChIP_Seq_and_RNA_Seq_overlap$entrez
head(CR.foldchanges)

keggres.sigmet.idx = gage(CR.foldchanges, gsets=kegg.sigmet.idx, same.dir=TRUE)
lapply(keggres.sigmet.idx, head,10)
fc.go.bp.p <- gage(CR.foldchanges, gsets = go.bp.gs,same.dir=TRUE)
lapply(fc.go.bp.p, head,20)
write.table(fc.go.bp.p$greater, file = "ChIP_Seq_and_RNA_Seq_overlap.go.bp.p.greater.txt",sep = "\t")
write.table(fc.go.bp.p$less, file = "ChIP_Seq_and_RNA_Seq_overlap.go.bp.p.less.txt",sep = "\t")
}
