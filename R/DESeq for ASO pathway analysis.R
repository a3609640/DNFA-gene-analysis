####################
## 1 Preparations ##
####################
# set global chunk options and load the neccessary packages
chooseCRANmirror()

source("http://bioconductor.org/biocLite.R")
biocLite("genefilter")
biocLite("apeglm")
biocLite("BiocStyle")
biocLite("rmarkdown")
biocLite("DESeq2")
biocLite("pathview")
biocLite("gage")
biocLite("gageData")
biocLite("RcppArmadillo")
biocLite("ReactomePA")
install.packages("https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.6.100.0.0.tar.gz", repos=NULL, type="source")

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


################################
## 2 Preparing count matrices ##
################################
# obtain the count table of the experiment directly from a pre-saved file. The RNA-seq was aligned to Hg38 by STAR
# read RNA-seq read data. 
setwd("~/Documents/Su Wu/Documents/Research/Naar Lab/RNA-seq/Test run/Analysis/STAR/results")
testseq <- read.csv("testseq.csv")
# Use the column one (Ensemble names) as columnn names. 
testseq <- data.frame(testseq[,-1], row.names=testseq[,1])
# Remove the first four rows (N_unmapped,N_multimapping,N_noFeature and N_ambiguous)
testseq <- data.frame(testseq[c(-1,-2,-3,-4),])
par(mar=c(3,12,2,1))
boxplot(testseq, outline=FALSE, horizontal=TRUE, las=1)

###################################
###################################
## 3. Compare ASO-Neg and ASO-4 ###
###################################
###################################
## generate dataset for ASO4 treatment
<<<<<<< HEAD
testASO <- data.frame(testseq[,c(5,6,7,8,21,22,23,24)])
=======
testASO <- data.frame(testseq[,c(13,14,15,16,21,22,23,24)])
>>>>>>> df2833f82e6dd4c02d8f5230e1d1d7b59bfb706c

par(mar=c(3,12,2,1))
boxplot(testASO, outline=FALSE, horizontal=TRUE, las=1)

## Prefiltering: by removing rows in which there are no reads or nearly no reads
guideDataASO <- testASO[rowSums(testASO)>1,]
head(guideDataASO)
dim(guideDataASO)
## Lets see how the data looks in a box plot
par(mar=c(3,12,2,1))
boxplot(guideDataASO, outline=FALSE, horizontal=TRUE, las=1)

# Time to create a design for our "modelling" 
guideDesignASO <- data.frame(row.names = colnames(guideDataASO),
<<<<<<< HEAD
                             condition = c(rep("siNeg", 4),
                                         #  rep("siSREBF1", 4),
                                           rep("ASO-4", 4)))
=======
                             condition = c(rep("ASO-Neg",4),rep("ASO-4",4)))
>>>>>>> df2833f82e6dd4c02d8f5230e1d1d7b59bfb706c

## object construction
## Construct DESeqDataSet with the count matrix, countData, and the sample information, colData
ddsASO <- DESeqDataSetFromMatrix(countData = guideDataASO,colData = guideDesignASO,design = ~ condition)
ddsASO

########################################
## 4 Exploring and exporting results ##
########################################


###############################################################################
## 4.1 standard analysis to make MA-plot from base means and log fold changes##
###############################################################################
<<<<<<< HEAD
ddsDEASO <- DESeq(ddsASO, test = "LRT", reduced = ~1)
=======
ddsDEASO <- DESeq(ddsASO)
>>>>>>> df2833f82e6dd4c02d8f5230e1d1d7b59bfb706c
resASO<-results(ddsDEASO) # default alpha = 0.1
## p-values and adjusted p-values
## We can order our results table by the smallest p value:
# ressiRNA<-ressiRNA[order(ressiRNA$log2FoldChange),]
resASO <- resASO[order(resASO$pvalue),]
## We can summarize some basic tallies using the summary function
summary(resASO)
## How many adjusted p-values were less than 0.1?
sum(resASO$padj < 0.1, na.rm=TRUE)

## Note that the results function automatically performs independent filtering based on the mean of normalized counts for each gene, 
## optimizing the number of genes which will have an adjusted p value below a given FDR cutoff, alpha = 0.05 here
res05 <- results(ddsDEASO, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

## plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. 
## Points will be colored red if the adjusted p value is less than 0.1. 
## Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(resASO, main="DESeq2 ASO", ylim=c(-2,2))
## It is more useful visualize the MA-plot for the shrunken log2 fold changes, 
## which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
resultsNames(ddsDEASO)
# normal is the the original DESeq2 shrinkage estimator, an adaptive normal prior
resLFC <- lfcShrink(resASO, coef=2)
# par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="normal")
# Alternative shrinkage estimators
# apeglm is the adaptive t prior shrinkage estimator from the apeglm package
resApe <- lfcShrink(resASO, coef=2, type="apeglm")
plotMA(resApe, xlim=xlim, ylim=ylim, main="apeglm")
# ashr is the adaptive shrinkage estimator from the ashr package (Stephens 2016). 
# Here DESeq2 uses the ashr option to fit a mixture of normal distributions to form the prior, with method="shrinkage"
resAsh <- lfcShrink(resASO, coef=2, type="ashr")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")


###########################################################
## 4.2 Add Entrez IDs, gene symbols, and full gene names ##
###########################################################
resASO <- data.frame(resASO)
columns(org.Hs.eg.db)
resASO$symbol = mapIds(org.Hs.eg.db,
                         keys=row.names(resASO), 
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
resASO$entrez = mapIds(org.Hs.eg.db,
                         keys=row.names(resASO), 
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
resASO$name =   mapIds(org.Hs.eg.db,
                         keys=row.names(resASO), 
                         column="GENENAME",
                         keytype="ENSEMBL",
                         multiVals="first")
summary(resASO)
head(resASO, 10)

########################################
## 4.3 Exporting results to CSV files ##
########################################
## Exporting only the results which pass an adjusted p value threshold can be accomplished with the subset function, 
## followed by the write.csv function.
resSig <- subset(resASO, padj < 0.1)
resSig
<<<<<<< HEAD
#write.csv(as.data.frame(resSig), 
#          file="ASO4_Neg_sigresults.csv")
# write.csv(as.data.frame(ressiRNA), 
#          file="ASO4_Neg_results.csv")
=======
write.csv(as.data.frame(resSig), 
          file="ASO4_Neg_sigresults.csv")
write.csv(as.data.frame(ressiRNA), 
          file="ASO4_Neg_results.csv")
>>>>>>> df2833f82e6dd4c02d8f5230e1d1d7b59bfb706c

#################################
## 4.4 KEGG pathways analysis ###
#################################
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
foldchanges = -resASO$log2FoldChange
names(foldchanges) = resASO$entrez
head(foldchanges)

# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
keggres.sigmet = gage(foldchanges, gsets=kegg.sets.hs.sigmet, same.dir=TRUE)
gores = gage(foldchanges, gsets=go.gs, same.dir=TRUE)
cartares = gage(foldchanges, gsets=carta.gs, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head,10)
lapply(keggres.sigmet, head,10)
lapply(gores, head,10)
lapply(cartares, head,10)

# Get KEGG pathway with only metabolism and signaling pathways
kg.hsa=kegg.gsets("hsa")
kegg.sigmet.idx=kg.hsa$kg.sets[kg.hsa$sigmet.idx]
keggres.sigmet.idx = gage(foldchanges, gsets=kegg.sigmet.idx, same.dir=TRUE)
<<<<<<< HEAD
lapply(keggres.sigmet.idx, head,10)
=======
lapply(keggres.sigmet.idx, head,20)
>>>>>>> df2833f82e6dd4c02d8f5230e1d1d7b59bfb706c
# write.table(keggres.sigmet.idx$greater, file = "keggres.sigmet.idx.greater.txt",sep = "\t")
# write.table(keggres.sigmet.idx$less, file = "keggres.sigmet.idx.less.txt",sep = "\t")

kegg.met.idx=kg.hsa$kg.sets[kg.hsa$met.idx]
keggres.met.idx = gage(foldchanges, gsets=kegg.met.idx, same.dir=TRUE)
lapply(keggres.met.idx, head,10)

kegg.sig.idx=kg.hsa$kg.sets[kg.hsa$sig.idx]
keggres.sig.idx = gage(foldchanges, gsets=kegg.sig.idx, same.dir=TRUE)
lapply(keggres.sig.idx, head,10)
