#####################
## 1. Preparations ##
#####################
# set global chunk options and load the neccessary packages
source("http://bioconductor.org/biocLite.R")
biocLite("genefilter")
biocLite("BiocStyle")
biocLite("rmarkdown")
biocLite("DESeq2")
biocLite("pathview")
biocLite("gage")
biocLite("gageData")
biocLite("RcppArmadillo")
biocLite("plotly")
biocLite("FactoMineR")
biocLite("factoextra")
install.packages("https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.6.100.0.0.tar.gz", 
                 repos=NULL, type="source")

library(FactoMineR)
library(factoextra)
library(cluster)
library(pheatmap)
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
library(SummarizedExperiment)
library(plotly)

#################################
## 2. Preparing count matrices ##
#################################
# obtain the count table of the experiment directly from a pre-saved file: testseq.csv. 
# The RNA-seq was aligned to human reference genome Hg38 by STAR aligner
# read processed RNA-seq read data from file testseq.csv. 
setwd("~/Documents/Su Wu/Documents/Research/Naar Lab/RNA-seq/Test run/Analysis/STAR/results")
testseq <- read.csv("testseq.csv")
# Use the column one (Ensemble names) as columnn names. 
testseq <- data.frame(testseq[,-1], row.names=testseq[,1])
# Remove the first four rows (N_unmapped,N_multimapping,N_noFeature and N_ambiguous)
testseq <- data.frame(testseq[c(-1,-2,-3,-4),])

##########################################
## 3. Quality control of the count data ##
##########################################
## check the read distribution of RNA-Seq results
par(mar=c(3,12,2,1))
boxplot(testseq, outline=FALSE, horizontal=TRUE, las=1)
## Remove rows in which there are no reads or nearly no reads
guideData <- testseq[rowSums(testseq)>1,]
head(guideData)
dim(guideData)
## check how the data looks in a box plot after removing rows with no read
par(mar=c(3,12,2,1))
boxplot(guideData, outline=FALSE, horizontal=TRUE, las=1)
## create a design for our "modelling" 
## each sample contains four techinical replicates
condition = c(rep("Mock",4),rep("siNegative",4),rep("siSREBF1",4),
              rep("ASO-Neg",4),rep("ASO-1",4),rep("ASO-4",4))
guideDesign <- data.frame(row.names = colnames(guideData),
                          condition = c(rep("Mock",4),rep("siNegative",4),rep("siSREBF1",4),
                                        rep("ASO-Neg",4),rep("ASO-1",4),rep("ASO-4",4)))
## Construct DESeqDataSet with the count matrix, countData, and the sample information, colData
dds <- DESeqDataSetFromMatrix(countData = guideData,colData = guideDesign,design = ~ condition)
dds
head(assay(dds))

##########################
## 4. Standard analysis ##
##########################
# DESeq function performs a default analysis through the steps:
# (1) estimation of size factor: estimateSizeFactors
# (2) estimation of dispersion: estimateDispersions
# (3) Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
ddsDE <- DESeq(dds)  
ddsres <- results(ddsDE)  
summary(ddsres)
res <- data.frame(ddsres)

######################################################
## 5. Regularized log transformation for clustering ##
######################################################
## The regularized log transform can be obtained using the rlog() function. 
## Regularized log transform is to stabilize the variance of the data and to make its distribution roughly symmetric
## The default “blinds” the normalization to the design. 
## This is very important so as to not bias the analyses (e.g. class discovery)
## The running times are shorter when using blind=FALSE and if the function DESeq has already been run, 
## because then it is not necessary to re-estimate the dispersion values. 
## The assay function is used to extract the matrix of normalized value
vsd <- vst(ddsDE, blind=FALSE)
rld <- rlog(ddsDE,blind=FALSE)
# Hierarchical clustering using rlog transformation
dists=dist(t(assay(rld)))
plot(hclust(dists), labels=guideDesign$condition)
sampleDistMatrix <- as.matrix(dists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=dists,
         clustering_distance_cols=dists,
         col=colors)

#####################################
## 6. Principal Component analysis ##
#####################################
## number of top genes to use for principal components, selected by highest row variance, 500 by default
data <- plotPCA(rld, intgroup = c( "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
## Print 2D PCA plot
black.bold.18.text <- element_text(face = "bold", color = "black", size = 18)
ggplot(data=data, aes_string(x="PC1", y="PC2", color="condition")) + 
      geom_point(size=3) + 
      theme_bw() + 
      xlim(-10, 6) + 
      ylim(-6, 6) +
      theme(text = black.bold.18.text, 
            axis.text = black.bold.18.text,
            axis.line.x = element_line(color="black", size=1),
            axis.line.y = element_line(color="black", size=1),
            axis.ticks = element_line(size = 1),
            axis.ticks.length = unit(.25, "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black",size=1),
            panel.background = element_blank(),
            legend.position=c(0,0),
            legend.justification=c(-0.05,-0.05)) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) 

## Print 3D PCA plot
plotPCA3D <- function (object, intgroup = "condition", ntop = 5000, returnData = FALSE){
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1],
                  PC2 = pca$x[, 2],
                  PC3 = pca$x[, 3],
                  group = group,
                  intgroup.df,
                  name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:3]
    return(d)
  }
  message("Generating plotly plot")
  p <- plotly::plot_ly(data = d,
                       x = ~PC1,
                       y = ~PC2,
                       z = ~PC3,
                       color = group,
                       mode = "markers",
                       type = "scatter3d")
  return(p)
}

plotPCA3D(rld, intgroup = "condition", ntop = 5000, returnData = FALSE)

##########################################################
## 7. Add Entrez IDs, gene symbols, and full gene names ##
##########################################################
columns(org.Hs.eg.db)
res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")
summary(res)
head(res, 10)


######################################################
## 8. Find the genes enriched in each PCA component ##
######################################################
## scale.unit : a logical value. If TRUE, the data are scaled to unit variance before the analysis. 
## This standardization to the same scale avoids some variables to become dominant just because of their large measurement units.
## We used FAlSE for scale.unit because rld has been run with DESEQ function before. 
## ncp : number of dimensions kept in the final results.
## graph : a logical value. If TRUE a graph is displayed.
head(assay(rld))
assayrld <- assay(rld)
Pvars <- rowVars(assayrld)
select <- order(Pvars, decreasing = TRUE)[seq_len(min(500, 
                                                      length(Pvars)))]

columns(org.Hs.eg.db)
row.names(assayrld) = mapIds(org.Hs.eg.db,
                    keys=row.names(assayrld), 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")


assayrld <-data.frame(t(assayrld[select, ]))
assayrld$condition = guideDesign$condition

con = c("Mock-1", "Mock-2", "Mock-3", "Mock-4",
        "siNegative-1", "siNegative-2", "siNegative-3", "siNegative-4",
        "siSREBF1-1", "siSREBF1-2", "siSREBF1-3", "siSREBF1-4",
        "ASO-Neg-1", "ASO-Neg-2", "ASO-Neg-3", "ASO-Neg-4",
        "ASO-1-1", "ASO-1-2","ASO-1-3","ASO-1-4",
        "ASO-4-1", "ASO-4-2", "ASO-4-3", "ASO-4-4")
rownames(assayrld) = con

# The variable Species (index = 501) is removed
# before PCA analysis
## # Compute PCA with ncp = 3, to keep only the first three principal components
res.pca <- PCA(assayrld[,-501], scale.unit = FALSE, ncp = 2,graph = TRUE)

#################################################################################
## 8. Extract the proportion of variances retained by the principal components ## 
#################################################################################
## Eigenvalues correspond to the amount of the variation explained by each principal component (PC). 
## Eigenvalues are large for the first PC and small for the subsequent PCs.
eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])
eigen <- eigenvalues[1:10,]
# Make a scree plot using base graphics : 
# A scree plot is a graph of the eigenvalues/variances associated with components.
barplot(eigen[, 2], names.arg=1:nrow(eigen), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")

lines(x = 1:nrow(eigen), eigen[, 2], 
      type="b", pch=19, col = "red")

# plot biplot PCA graph with the top six contributors for PC1 and PC2
fviz_pca_biplot(res.pca, 
                select.var = list(contrib = 6),
                #select.var = list(contrib = 0.6), 
                col.var = "red",
                label="var", 
                habillage=assayrld$condition)+ 
  geom_point(size=3, 
             aes(colour = factor(assayrld$condition))) + 
  theme_bw() + 
  xlim(-8, 4) + 
  ylim(-5, 5) +
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank(),
        legend.text = element_text(colour="black", size = 18, face = "bold"),  
        legend.position=c(0,1),
        legend.justification=c(-0.05,1.05)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) 

########################################################
## 9. Hierarchical Clustering on Principal Components ## 
########################################################
## Compute hierarchical clustering: Hierarchical clustering is performed using the Ward’s criterion on the selected principal components. 
## Ward criterion is used in the hierarchical clustering because it is based on the multidimensional variance like principal component analysis.
res.hcpc <- HCPC(res.pca, graph = FALSE)
fviz_dend(res.hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
          )
fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
            )


#Hierarchical Clustering

#compute dissimilarity matrix for all data
eu.d <- dist(data, method = "euclidean")
# Hierarchical clustering using Ward's method
res.hc <- hclust(eu.d, method = "ward.D2" )
# Cut tree into 4 groups
grp <- cutree(res.hc, k = 2)
# Visualize
plot(res.hc, cex = 0.6) # plot tree
rect.hclust(res.hc, k = 2, border = c("yellow","blue")) # add rectangle

## Contributions of variables to PCs
head(res.pca$var$contrib,10)
head(res.pca$var$cos2, 10)
var <- get_pca_var(res.pca)
var
library("corrplot")
corrplot(var$contrib, is.corr=FALSE)   
# Contributions of variables to PC1
fviz_contrib(res.pca, choice="var", axes = 1,top = 10, 
             fill = "lightgray", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45))
# Contributions of variables to PC2
fviz_contrib(res.pca, choice="var", axes = 2,top = 10, 
             fill = "lightgray", color = "black") +
             theme_minimal() +
             theme(axis.text.x = element_text(angle=45))

##  identify the most correlated variables with a given principal component
res.desc <- dimdesc(res.pca, axes = c(1,2),proba = 0.05)
# Description of dimension 1
head(res.desc, 10)
res.desc$Dim.1
dim1 <-data.frame(res.desc$Dim.1)
columns(org.Hs.eg.db)
dim1$symbol = mapIds(org.Hs.eg.db,
                     keys=row.names(dim1), 
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
dim1$entrez = mapIds(org.Hs.eg.db,
                     keys=row.names(dim1), 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
dim1$name =   mapIds(org.Hs.eg.db,
                     keys=row.names(dim1), 
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")
summary(dim1)
head(dim1,10)

res.desc$Dim.2
dim2 <-data.frame(res.desc$Dim.2)
columns(org.Hs.eg.db)
dim2$symbol = mapIds(org.Hs.eg.db,
                     keys=row.names(dim2), 
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
dim2$entrez = mapIds(org.Hs.eg.db,
                     keys=row.names(dim2), 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
dim2$name =   mapIds(org.Hs.eg.db,
                     keys=row.names(dim2), 
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")
summary(dim2)
head(dim2,10)

quanti.correlation.dim1 = dim1$quanti.correlation 
names(quanti.correlation.dim1) = dim1$entrez
head(quanti.correlation.dim1)

keggres.1 = gage(quanti.correlation.dim1, gsets=kegg.sets.hs, same.dir=TRUE)
keggres.sigmet.1 = gage(quanti.correlation.dim1, gsets=kegg.sets.hs.sigmet, same.dir=TRUE)
gores.1 = gage(quanti.correlation.dim1, gsets=go.gs, same.dir=TRUE)
cartares.1 = gage(quanti.correlation.dim1, gsets=carta.gs, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
lapply(keggres.1, head,10)
lapply(keggres.sigmet.1, head,10)
lapply(gores.1, head,10)
lapply(cartares.1, head,10)



quanti.correlation.dim2 = dim2$quanti.correlation 
names(quanti.correlation.dim2) = dim2$entrez
head(quanti.correlation.dim2)

keggres.2 = gage(quanti.correlation.dim2, gsets=kegg.sets.hs, same.dir=TRUE)
keggres.sigmet.2 = gage(quanti.correlation.dim2, gsets=kegg.sets.hs.sigmet, same.dir=TRUE)
gores.2 = gage(quanti.correlation.dim2, gsets=go.gs, same.dir=TRUE)
cartares.2 = gage(quanti.correlation.dim2, gsets=carta.gs, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
lapply(keggres.2, head,10)
lapply(keggres.sigmet.2, head,10)
lapply(gores.2, head,10)
lapply(cartares.2, head,10)


##################################################################
## 10. Plot of normalized counts for a single gene on log scale ##
##################################################################
# plotcount: "normalized" whether the counts should be normalized by size factor (default is TRUE)
# plotcount: "transform" whether to present log2 counts (TRUE) or to present the counts on the log scale (FALSE, default)
# re-arrange x-ase according to the following order: "Mock","siNeg","siBF1","ASO-neg","ASO-1","ASO-4"
# theme_bw() removes background color in the graph, guides(fill=FALSE) removes legends
black.bold.18.text <- element_text(face = "bold", color = "black", size = 18)

SREBF1 <- plotCounts(dds, gene="ENSG00000072310", 
                     intgroup="condition", 
                     normalized = TRUE, 
                     transform = TRUE,
                     returnData=TRUE)
SREBF1$condition <- factor(SREBF1$condition, 
                           levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(SREBF1, aes(x=condition, y=log2(count), fill=condition)) + 
      geom_boxplot()+ 
      ylim(8, 10.5)+ 
      theme_bw()+ 
      guides(fill=FALSE)+
      theme(text = black.bold.18.text, 
            axis.text = black.bold.18.text,
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.line.x = element_line(color="black", size=1),
            axis.line.y = element_line(color="black", size=1),
            axis.ticks = element_line(size = 1),
            axis.ticks.length = unit(.25, "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black",size=1),
            panel.background = element_blank())+
      labs(title = "SREBF1",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
SCD <- plotCounts(dds, gene="ENSG00000099194", 
                     intgroup="condition", 
                     normalized = TRUE, 
                     transform = TRUE,
                     returnData=TRUE)
SCD$condition <- factor(SCD$condition, 
                           levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(SCD, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(13, 14.5)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "SCD",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
FASN <- plotCounts(dds, gene="ENSG00000169710", 
                  intgroup="condition", 
                  normalized = TRUE, 
                  transform = TRUE,
                  returnData=TRUE)
FASN$condition <- factor(FASN$condition, 
                        levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(FASN, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(12, 13.25)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "FASN",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
ACACA <- plotCounts(dds, gene="ENSG00000278540", 
                  intgroup="condition", 
                  normalized = TRUE, 
                  transform = TRUE,
                  returnData=TRUE)
ACACA$condition <- factor(ACACA$condition, 
                        levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(ACACA, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(9.7, 10.5)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "ACACA",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
ACSL1 <- plotCounts(dds, gene="ENSG00000169710", 
                   intgroup="condition", 
                   normalized = TRUE, 
                   transform = TRUE,
                   returnData=TRUE)
ACSL1$condition <- factor(ACSL1$condition, 
                         levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(ACSL1, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(11.5, 13.5)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "ACSL1",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
BRAF <- plotCounts(dds, gene="ENSG00000157764", 
                  intgroup="condition", 
                  normalized = TRUE, 
                  transform = TRUE,
                  returnData=TRUE)
BRAF$condition <- factor(SCD$condition, 
                        levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(BRAF, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(6, 9)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "BRAF",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
SPIRE1 <- plotCounts(dds, gene="ENSG00000134278", 
                   intgroup="condition", 
                   normalized = TRUE, 
                   transform = TRUE,
                   returnData=TRUE)
SPIRE1$condition <- factor(SPIRE1$condition, 
                         levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(SPIRE1, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(4, 9)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "SPIRE1",x=" ", y= "log2(read counts)")
## *************************************************************
ENSG00000124486
## *************************************************************
USP9X <- plotCounts(dds, gene="ENSG00000124486", 
                     intgroup="condition", 
                     normalized = TRUE, 
                     transform = TRUE,
                     returnData=TRUE)
USP9X$condition <- factor(USP9X$condition, 
                           levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(USP9X, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(7, 10)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "USP9X",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
INSIG1 <- plotCounts(dds, gene="ENSG00000186480", 
                   intgroup="condition", 
                   normalized = TRUE, 
                   transform = TRUE,
                   returnData=TRUE)
INSIG1$condition <- factor(INSIG1$condition, 
                         levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(INSIG1, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(10, 12)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "INSIG1",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
YY1 <- plotCounts(dds, gene="ENSG00000100811", 
                     intgroup="condition", 
                     normalized = TRUE, 
                     transform = TRUE,
                     returnData=TRUE)
YY1$condition <- factor(YY1$condition, 
                           levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(YY1, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(8, 10.5)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "YY1",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
FABP7 <- plotCounts(dds, gene="ENSG00000164434", 
                  intgroup="condition", 
                  normalized = TRUE, 
                  transform = TRUE,
                  returnData=TRUE)
FABP7$condition <- factor(FABP7$condition, 
                        levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(FABP7, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(6, 9)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "FABP7",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
MITF <- plotCounts(dds, gene="ENSG00000187098", 
                   intgroup="condition", 
                   normalized = TRUE, 
                   transform = TRUE,
                   returnData=TRUE)
MITF$condition <- factor(MITF$condition, 
                         levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(MITF, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(10, 12)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "MITF",x=" ", y= "log2(read counts)")

## *************************************************************
SREBF2 <- plotCounts(dds, gene="ENSG00000198911", 
                   intgroup="condition", 
                   normalized = TRUE, 
                   transform = TRUE,
                   returnData=TRUE)
SREBF2$condition <- factor(SREBF2$condition, 
                         levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(SREBF2, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(10, 12)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "SREBF2",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
HMGCR <- plotCounts(dds, gene="ENSG00000113161", 
                     intgroup="condition", 
                     normalized = TRUE, 
                     transform = TRUE,
                     returnData=TRUE)
HMGCR$condition <- factor(HMGCR$condition, 
                           levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(HMGCR, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(10, 12)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "HMGCR",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
NRAS <- plotCounts(dds, gene="ENSG00000213281", 
                    intgroup="condition", 
                    normalized = TRUE, 
                    transform = TRUE,
                    returnData=TRUE)
NRAS$condition <- factor(NRAS$condition, 
                          levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(NRAS, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(7, 9)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "NRAS",x=" ", y= "log2(read counts)")
## *************************************************************




#############################
#############################
## Compare siNeg and siBF1 ##
#############################
#############################
## generate dataset for siRNA treatment
testsiRNA <- data.frame(testseq[,c(5,6,7,8,9,10,11,12)])

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
ddssiRNA <- DESeqDataSetFromMatrix(countData = guideDatasiRNA,colData = guideDesignsiRNA,design = ~ condition)
ddssiRNA

###########################################################################
## standard analysis to make MA-plot from base means and log fold changes##
###########################################################################
ddsDEsiRNA <- DESeq(ddssiRNA)
ressiRNA<-results(ddsDEsiRNA)
ressiRNA<-ressiRNA[order(ressiRNA$log2FoldChange),]
## An MA-plot21 provides a useful overview for an experiment with a two-group comparison
## The log2 fold change for a particular comparison is plotted on the y-axis and the average of the counts normalized by size factor is shown on the x-axis 
## (“M” for minus, because a log ratio is equal to log minus log, and “A” for average). Each gene is represented with a dot. 
## Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red.
plotMA(ddsDEsiRNA, main="DESeq2", ylim=c(-2,2))
## The red points indicate genes for which the log2 fold change was significantly higher than 0.5 or less than -0.5 
# (treatment resulting in more than doubling or less than halving of the normalized counts) with adjusted p value less than 0.1.
resLFC1 <- results(ddsDEsiRNA, lfcThreshold=0.5)
table(resLFC1$padj < 0.1)
plotMA(resLFC1,main="DESeq2", ylim=c(-2,2))

#######################################################
## Add Entrez IDs, gene symbols, and full gene names ##
#######################################################
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
topVarGenes <- head(order(rowVars(assay(vdsiRNA)), decreasing=TRUE), 200)
mat <- assay(rldsiRNA)[topVarGenes, ]
mat <- mat - rowMeans(mat)
heatmap.2(mat,labRow = NA,labCol=NA, scale="row",lhei = c(2, 8),
          trace="none", dendrogram="none", margins=c(2,10),
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          ColSideColors = c(Control="gray", DPN="darkgreen", OHT="orange")[colData(rldsiRNA)$condition ])


############################
## KEGG pathways analysis ##
############################

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
foldchanges = -ressiRNA$log2FoldChange
names(foldchanges) = ressiRNA$entrez
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

# Get the pathways
library(dplyr)
keggrespathways <- data.frame(id=rownames(keggres$less), keggres$less) %>%
  tbl_df() %>%
  filter(row_number()<=10) %>%
  .$id %>%
  as.character()

keggrespathways

# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids


# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
# plot multiple pathways (plots saved to disk and returns a throwaway list object)
detach("package:dplyr",unload=TRUE)
pv.out.list <- sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))

write.table(keggres.sigmet$greater, file = "keggres.sigmet.greater.txt",sep = "\t")
write.table(keggres.sigmet$less, file = "keggres.sigmet.less.txt",sep = "\t")

## map metabolites on a large scale using pathway.id= '01100'
## Cancer  - Homo sapiens (human) pathway.id= '05200'
## Melanoma - Homo sapiens (human) pathway.id= '05218'
## hsa04010 MAPK signaling pathway   
## hsa04910 Insulin signaling pathway 
## hsa00061 fatty acid biosynthesis pathway 
pv.out <- pathview(gene.data = foldchanges, pathway.id= 'hsa04210', species = "hsa", kegg.native = T, same.layer = F)


######################################################
## derive the non-redundant signficant gene set lists
kegg.p <- gage(foldchanges, gsets=kegg.sets.hs.sigmet,same.dir=TRUE)
kegg.esg.up <- esset.grp(kegg.p$greater,
                         foldchanges, gsets=kegg.sets.hs.sigmet,
                         test4up = TRUE, output = TRUE, outname = "kegg.up", make.plot = FALSE)

names(kegg.esg.up)
head(kegg.p,4)
head(kegg.esg.up$essentialSets, 4)
head(kegg.esg.up$setGroups, 4)
head(kegg.esg.up$coreGeneSets, 40)
######################################################
## kegg test for 2-directional changes
kegg.2d.p <- gage(foldchanges, gsets = kegg.sets.hs.sigmet,same.dir = FALSE)
head(kegg.2d.p$greater)
head(kegg.2d.p$stats)
rownames(kegg.2d.p$greater)[1:3]

#################################
## Off-target analysis of ASO  ##
#################################
## generate dataset for siRNA treatment
testASO <- data.frame(testseq[,c(9,10,11,12,17,18,19,20)])

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
                             condition = c(rep("siBF1",4),rep("ASO-BF1",4)))

## object construction
## Construct DESeqDataSet with the count matrix, countData, and the sample information, colData
ddsASO <- DESeqDataSetFromMatrix(countData = guideDataASO,colData = guideDesignASO,design = ~ condition)
ddsASO

###########################################################################
## standard analysis to make MA-plot from base means and log fold changes##
###########################################################################
ddsDEASO <- DESeq(ddsASO)
plotMA(ddsDEASO, main="DESeq2", ylim=c(-2,2))
ddsresASO <- results(ddsDEASO)
plotMA(ddsresASO,main="DESeq2", ylim=c(-2,2))

# results extracts a result table from a DESeq analysis giving base means across samples, 
# log2 fold changes, standard errors, test statistics, p-values and adjusted p-values.
# lfcThreshold:a non-negative value which specifies a log2 fold change threshold. 
# The default value is 0, corresponding to a test that the log2 fold changes are equal to zero. 
# The user can specify the alternative hypothesis using the altHypothesis argument, which defaults to testing for log2 fold changes greater in absolute value than a given threshold. 

resLFC1ASO <- results(ddsDEASO, lfcThreshold=1)
table(resLFC1ASO$padj < 0.1)
plotMA(resLFC1ASO,main="DESeq2", ylim=c(-2,2))

# a scatter plot of log2 fold changes (on the y-axis) versus the mean of normalized counts (on the x-axis)
plotMA(resLFC1ASO, main="siSREBF1 vs ASO-1",ylim=c(-3,3), xlab="Mean of Normalized Counts", ylab="Log2 Fold Change")
topGene <- rownames(resLFC1ASO)[which.max(resLFC1ASO$log2FoldChange)]
with(resLFC1ASO[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})



resASO = ddsresASO[order(ddsresASO$pvalue),]
resASO <- data.frame(resASO)
summary(resASO)

###########################################################
## Regularized log transformation for clustering and PCA ##
###########################################################
# The regularized log transform can be obtained using the rlog() function. 
# Note that an important argument for this function is blind (TRUE by default). 
# The default “blinds” the normalization to the design. 
# This is very important so as to not bias the analyses (e.g. class discovery)
rldASO=rlog(ddsASO,blind=TRUE)
topVarGenes <- head( order( rowVars( assay(rldASO) ), decreasing=TRUE ), 50 )
heatmap.2( assay(rldASO)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", margins=c(2,10),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" )[
             colData(rldASO)$condition ] )

#######################################################
## Add Entrez IDs, gene symbols, and full gene names ##
#######################################################
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

#############################
#############################
## Compare ASO-Neg and ASO-4 ##
#############################
#############################
## generate dataset for siRNA treatment
testASO <- data.frame(testseq[,c(13,14,15,16,21,22,23,24)])

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
                               condition = c(rep("ASO-Neg",4),rep("ASO-4",4)))

## object construction
## Construct DESeqDataSet with the count matrix, countData, and the sample information, colData
ddsASO <- DESeqDataSetFromMatrix(countData = guideDataASO,colData = guideDesignASO,design = ~ condition)
ddsASO

###########################################################################
## standard analysis to make MA-plot from base means and log fold changes##
###########################################################################
ddsDEASO <- DESeq(ddsASO)
resASO<-results(ddsDEASO)
resASO<-resASO[order(resASO$log2FoldChange),]
## An MA-plot21 provides a useful overview for an experiment with a two-group comparison
## The log2 fold change for a particular comparison is plotted on the y-axis and the average of the counts normalized by size factor is shown on the x-axis 
## (“M” for minus, because a log ratio is equal to log minus log, and “A” for average). Each gene is represented with a dot. 
## Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red.
plotMA(ddsDEASO, main="DESeq2 ASO", ylim=c(-2,2))
## The red points indicate genes for which the log2 fold change was significantly higher than 0.5 or less than -0.5 
# (treatment resulting in more than doubling or less than halving of the normalized counts) with adjusted p value less than 0.1.
resultsNames(ddsDEASO)
resLFC <- lfcShrink(ddsDEASO, coef="condition_ASO.Neg_vs_ASO.4", type="apeglm")
resLFC
plotMA(resLFC1,main="DESeq2", ylim=c(-2,2))

#######################################################
## Add Entrez IDs, gene symbols, and full gene names ##
#######################################################
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

############################
## KEGG pathways analysis ##
############################

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
kegg = gage(foldchanges, gsets=kegg.gs, same.dir=TRUE)
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
keggres.sigmet = gage(foldchanges, gsets=kegg.sets.hs.sigmet, same.dir=TRUE)
gores = gage(foldchanges, gsets=go.gs, same.dir=TRUE)
cartares = gage(foldchanges, gsets=carta.gs, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
lapply(kegg, head,10)
lapply(keggres, head,10)
lapply(keggres.sigmet, head,10)
lapply(gores, head,10)
lapply(cartares, head,10)


kg.hsa=kegg.gsets("hsa")
kegg.sigmet.idx=kg.hsa$kg.sets[kg.hsa$sigmet.idx]
keggres.sigmet.idx = gage(foldchanges, gsets=kegg.sigmet.idx, same.dir=TRUE)
lapply(keggres.sigmet.idx, head,10)

kegg.met.idx=kg.hsa$kg.sets[kg.hsa$met.idx]
keggres.met.idx = gage(foldchanges, gsets=kegg.met.idx, same.dir=TRUE)
lapply(keggres.met.idx, head,10)

kegg.sig.idx=kg.hsa$kg.sets[kg.hsa$sig.idx]
keggres.sig.idx = gage(foldchanges, gsets=kegg.sig.idx, same.dir=TRUE)
lapply(keggres.sig.idx, head,10)

