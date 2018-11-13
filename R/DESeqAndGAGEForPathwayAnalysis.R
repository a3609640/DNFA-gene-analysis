dataRoot <- file.path("project-data")

####################
## 1 Preparations ##
####################
# set global chunk options and load the neccessary packages
chooseCRANmirror()

#BiocManager::install("genefilter")
#BiocManager::install("apeglm")
#BiocManager::install("ashr")
#BiocManager::install("BiocStyle")
#BiocManager::install("biomaRt")
#BiocManager::install("DESeq2")
#BiocManager::install("gage")
#BiocManager::install("gageData")
#BiocManager::install("pathview")
#BiocManager::install("RcppArmadillo")
#BiocManager::install("ReactomePA")
#BiocManager::install("rmarkdown")
#install.packages("https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.6.100.0.0.tar.gz", repos=NULL, type="source")

library(ashr)
library(apeglm)
library(RcppArmadillo)
library(colorspace)
library(lattice)
library(RODBC)
library(Matrix)
library(survival)
library(Rcpp)
library(genefilter)
library(BiocStyle)
library(rmarkdown)
library(geneplotter)
library(ggplot2)
library(gplots)
library(plyr)
library(DESeq2)
library(RColorBrewer)
library(stringr)
library(biomaRt)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pathview)
library(gage)
library(gageData)
library(Biobase)
library(S4Vectors)
library(stats4)
library(BiocGenerics)
library(parallel)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(ReactomePA)
library(SummarizedExperiment)
#######

doAll2 <- function() {
################################
## 2 Preparing count matrices ##
################################
# obtain the count table of the experiment directly from a pre-saved file: gene-counts.csv. 
# The RNA-seq was aligned to human reference genome Hg38 by STAR aligner
# read processed RNA-seq read data from file testseq.csv.
testseqCSV <- file.path(dataRoot, "gene-counts-from-Makefile.csv")
testseq <- read.csv(testseqCSV)
# Use the column one (Ensemble names) as columnn names. 
testseq <- data.frame(testseq[,-1], row.names=testseq[,1])
# Remove the first four rows (N_unmapped,N_multimapping,N_noFeature and N_ambiguous)
testseq <- data.frame(testseq[c(-1,-2,-3,-4),])
par(mar=c(3,12,2,1))
boxplot(testseq, outline=FALSE, horizontal=TRUE, las=1)


###############################
## 3 Compare siNeg and siBF1 ##
###############################

## generate dataset for siRNA treatment
testsiRNA <- data.frame(testseq[,c(5,6,7,8,9,10,11,12)])
## check the read distribution by boxplot
par(mar=c(3,12,2,1))
boxplot(testsiRNA, outline=FALSE, horizontal=TRUE, las=1)

## Prefiltering: by removing rows in which there are no reads or nearly no reads
guideDatasiRNA <- testsiRNA[rowSums(testsiRNA)>1,]
head(guideDatasiRNA)
dim(guideDatasiRNA)
## Lets see how the data looks in a box plot
par(mar=c(3,12,2,1))
boxplot(guideDatasiRNA, outline=FALSE, horizontal=TRUE, las=1)

# Time to create a design for our "modelling" 
guideDesignsiRNA <- data.frame(row.names = colnames(guideDatasiRNA),
                               condition = c(rep("siNeg",4),rep("siBF1",4)))

## object construction
## Construct DESeqDataSet with the count matrix, countData, and the sample information, colData
ddssiRNA <- DESeqDataSetFromMatrix(countData = guideDatasiRNA,
                                   colData = guideDesignsiRNA,
                                   design = ~ condition)
## specifying the reference level: 
ddssiRNA$condition <- relevel(ddssiRNA$condition, ref = "siNeg")
ddssiRNA


########################################
## 4 Exploring and exporting results ##
########################################


###############################################################################
## 4.1 standard analysis to make MA-plot from base means and log fold changes##
###############################################################################
ddsDEsiRNA <- DESeq(ddssiRNA)
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


resApeT <- lfcShrink(ddsDEsiRNA, coef=2, lfcThreshold=.5)
plotMA(resApeT, ylim=c(-3,3), cex=.8)
abline(h=c(-0.5,0.5), col="dodgerblue", lwd=2)

# Alternative shrinkage estimators
# apeglm is the adaptive t prior shrinkage estimator from the apeglm package
resApe <- lfcShrink(ddsDEsiRNA, coef=2, type="apeglm")
plotMA(resApe, xlim=xlim, ylim=ylim, main="apeglm")
# ashr is the adaptive shrinkage estimator from the ashr package (Stephens 2016). 
# Here DESeq2 uses the ashr option to fit a mixture of normal distributions to form the prior, with method="shrinkage"
resAsh <- lfcShrink(ddsDEsiRNA, coef=2, type="ashr")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

###########################################################
## 4.2 Add Entrez IDs, gene symbols, and full gene names ##
###########################################################
resApe <- data.frame(resApe)
columns(org.Hs.eg.db)
resApe$symbol = mapIds(org.Hs.eg.db,
                       keys=row.names(resApe), 
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
resApe$entrez = mapIds(org.Hs.eg.db,
                       keys=row.names(resApe), 
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")
resApe$name =   mapIds(org.Hs.eg.db,
                       keys=row.names(resApe), 
                       column="GENENAME",
                       keytype="ENSEMBL",
                       multiVals="first")
summary(resApe)
head(resApe, 10)

###########################################################
###########################################################

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





#####################
### Volcano plot #### 
#####################
with(resApe, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(resApe, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(resApe, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(resApe, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(resApe, padj<.05 & abs(log2FoldChange)>1), 
     textxy(log2FoldChange, -log10(pvalue), labs=symbol, cex=.8))


head(ressiRNA)
with(ressiRNA, plot(log2FoldChange, -log10(pvalue), pch=20, main="eXpress DESeq2 Results", xlim=c(-5,5)))
axis(1, at = seq(-5, 5, by = 1))
with(subset(ressiRNA, padj<=0.05 & log2FoldChange<=-1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(ressiRNA, padj<=0.05 & log2FoldChange>=1), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
library(calibrate)
with(subset(ressiRNA, padj<.01 & abs(log2FoldChange)>2.32), textxy(log2FoldChange, -log10(pvalue), labs=symbol, cex=.5))




########################################
## 4.3 Exporting results to CSV files ##
########################################
## Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
## followed by the write.csv function.
resSig <- subset(ressiRNA, pvalue < 0.05)
resSig
write.csv(as.data.frame(resSig), 
          file="siBF1_siNeg_sigresults.csv")
write.csv(as.data.frame(ressiRNA), 
          file="siBF1_siNeg_results.csv")


###########################################################
## Heatmap for top variable genes across siNeg and siBF1 ##
###########################################################
# The regularized log transform can be obtained using the rlog() function. 
# The default “blinds” the normalization to the design. 
# topVarGenes looks at the row variance of the transformed values regardless of which samples come from which group. 
# Differential expression testing asks whether the difference across group is large relative to the within-group variance. 
# So these are different ways of ranking genes.
# calculate the variance for each gene, # select the top 50 genes by variance
rldsiRNA <- rlog(ddsDEsiRNA,blind=TRUE)
vdsiRNA <- varianceStabilizingTransformation(ddsDEsiRNA, blind=TRUE)
topVarGenes <- head(order(rowVars(assay(vdsiRNA)), decreasing=TRUE), 500)
mat <- assay(rldsiRNA)[topVarGenes, ]
mat <- mat - rowMeans(mat)
heatmap.2(mat,labRow = NA,labCol=NA, scale="row",lhei = c(2, 8),
          trace="none", dendrogram="none", margins=c(2,10),
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          ColSideColors = c(Control="gray", DPN="darkgreen", OHT="orange")[colData(rldsiRNA)$condition ])


###############################################
## KEGG pathways analysis on DESeq for siRNA ##
###############################################
## Generally Applicable Gene-set/Pathway Analysis (GAGE)
## check for coordinated differential expression over gene sets instead of changes of individual genes.
data(kegg.gs)
data(go.gs)
data(carta.gs)
lapply(kegg.gs[1:3],head)
lapply(go.gs[1:3],head)
## kegg.sets.hs is a named list of 229 elements
## Each element is a character vector of member gene Entrez IDs for a single KEGG pathway
data(kegg.sets.hs)
## sigmet.idx.hs is an index of numbers of sinaling and metabolic pathways in kegg.set.gs.
data(sigmet.idx.hs)
## kegg.sets.hs[sigmet.idx.hs] gives you the “cleaner” gene sets of sinaling and metabolic pathways only.
kegg.sets.hs.sigmet <-  kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs.sigmet, 3)

## Generate a named vector of fold changes, where the names of the values are the Entrez gene IDs, for the gage() function
foldchanges = ressiRNA$log2FoldChange
names(foldchanges) = ressiRNA$entrez
head(foldchanges)


# Get the results
# keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
# keggres.sigmet = gage(foldchanges, gsets=kegg.sets.hs.sigmet, same.dir=TRUE)
# gores = gage(foldchanges, gsets=go.gs, same.dir=TRUE)
# cartares = gage(foldchanges, gsets=carta.gs, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
# lapply(keggres, head,10)
# lapply(keggres.sigmet, head,10)
# lapply(gores, head,10)
# lapply(cartares, head,10)

# Get KEGG pathway with only metabolism and signaling pathways
kg.hsa=kegg.gsets("hsa")
kegg.sigmet.idx=kg.hsa$kg.sets[kg.hsa$sigmet.idx]
keggres.sigmet.idx = gage(foldchanges, gsets=kegg.sigmet.idx, same.dir=TRUE)
lapply(keggres.sigmet.idx, head,20)
write.table(keggres.sigmet.idx$greater, file = "keggres.sigmet.idx.greater.txt",sep = "\t")
write.table(keggres.sigmet.idx$less, file = "keggres.sigmet.idx.less.txt",sep = "\t")

kegg.met.idx=kg.hsa$kg.sets[kg.hsa$met.idx]
keggres.met.idx = gage(foldchanges, gsets=kegg.met.idx, same.dir=TRUE)
lapply(keggres.met.idx, head,10)

kegg.sig.idx=kg.hsa$kg.sets[kg.hsa$sig.idx]
keggres.sig.idx = gage(foldchanges, gsets=kegg.sig.idx, same.dir=TRUE)
lapply(keggres.sig.idx, head,10)

# set up GO database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP] # “Biological Process”
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF] # “Molecular Function”
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC] # “Cellular Component”

fc.go.bp.p <- gage(foldchanges, gsets = go.bp.gs,same.dir=TRUE)
lapply(fc.go.bp.p, head,20)
write.table(fc.go.bp.p$greater, file = "fc.go.bp.p.greater.txt",sep = "\t")
write.table(fc.go.bp.p$less, file = "fc.go.bp.p.less.txt",sep = "\t")

fc.go.mf.p <- gage(foldchanges, gsets = go.mf.gs)
lapply(fc.go.mf.p, head,10)
fc.go.cc.p <- gage(foldchanges, gsets = go.cc.gs)
lapply(fc.go.cc.p, head,10)

########################################################################
###### Pathway analysis on ChIP-seq and RNA-seq overlapping genes ######
########################################################################
setwd("~/Documents/Su Wu/Documents/Research/Naar Lab/ChIP-seq")
library(readxl)
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

