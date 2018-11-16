library(AnnotationDbi)
library(BiocStyle)
library(cluster)
library(colorspace)
library(DESeq2)
library(Biobase)
library(BiocGenerics)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(gage)
library(gageData)
library(geneplotter)
library(genefilter)
library(GenomeInfoDb)
library(GenomicRanges)
library(ggplot2)
library(gplots)
library(IRanges)
library(lattice)
library(Matrix)
library(org.Hs.eg.db)
library(parallel)
library(pheatmap)
library(plotly)
library(plyr)
library(rmarkdown)
library(Rcpp)
library(RcppArmadillo)
library(RColorBrewer)
library(S4Vectors)
library(SummarizedExperiment)
library(stats4)
library(stringr)
library(survival)

readGeneCounts <- function () {
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

##########################################################
## 2. Preparing count matrices from the RNA-seq results ##
##########################################################
getGuideData <- function() {
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

getGuideDesign <- function(guideData) {
  ## create a design for our "modelling" 
  ## each sample contains four techinical replicates
  #condition = c(rep("Mock",4),rep("siNegative",4),rep("siSREBF1",4),
  #              rep("ASO-Neg",4),rep("ASO-1",4),rep("ASO-4",4))
  return(data.frame(row.names = colnames(guideData),
                            condition = c(rep("Mock",4),rep("siNegative",4),rep("siSREBF1",4),
                                          rep("ASO-Neg",4),rep("ASO-1",4),rep("ASO-4",4))))
  
}

#######################################################
### 3. Construct DESeqDataSet from the count matrix ###
#######################################################
getDDS <- function (guideData) {
  ## Construct DESeqDataSet with the count matrix, countData, and the sample information, colData
  dds <- DESeqDataSetFromMatrix(countData = guideData,colData = guideDesign,design = ~ condition)
  dds
  head(assay(dds))
  return(dds)  
}

######################################################
#### 4. Standard differential expression analysis ####
######################################################
# DESeq function performs a default analysis through the steps:
# (1) estimation of size factor: estimateSizeFactors
# (2) estimation of dispersion: estimateDispersions
# (3) Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
getDDSRES <- function(ddsDE) {
  ddsres <- results(ddsDE)  
  summary(ddsres)
  res <- data.frame(ddsres)
  return(ddsres)  
}

########################################################
##### 5. Count data transformations for clustering #####
########################################################
## The regularized log transform can be obtained using the rlog() function. 
## Regularized log transform is to stabilize the variance of the data and to make its distribution roughly symmetric
## The default “blinds” the normalization to the design. 
## This is very important so as to not bias the analyses (e.g. class discovery)
## The running times are shorter when using blind=FALSE and if the function DESeq has already been run, 
## because then it is not necessary to re-estimate the dispersion values. 
## The assay function is used to extract the matrix of normalized value
makeHeatMap <- function(guideDesign, ddsDE) {
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
  return(rld)
}

doAll1 <- function() {
testseq <- readGeneCounts()
guideData <- getGuideData()
guideDesign <- getGuideDesign(guideData)
dds <- getDDS(guideData)
ddsDE <- DESeq(dds)
res <- getDDSRES(ddsDE)
rld <- makeHeatMap(guideDesign, ddsDE)

#############################################
###### 6. Principal component analysis ######
#############################################
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

####################################################################
####### 7. Add Entrez IDs, gene symbols, and full gene names #######
####################################################################
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

##################################################################
######## 8. Annotate genes enriched in each PCA component ########
##################################################################
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

####################################################################
######### 9. Extract variances in each principal component ######### 
####################################################################
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

# plot biplot graph with the top six contributing genes to PCA from RNA-Seq
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

#########################################################################
########## 10. Hierarchical Clustering on Principal Components ########## 
#########################################################################
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

##########################################################################
########### 11. Top contributing gene variables to PC1 and PC2 ###########
##########################################################################
## Contributions of variables to PCs
head(res.pca$var$contrib,10)
head(res.pca$var$cos2, 10)
var <- get_pca_var(res.pca)
var
#library("corrplot")
#corrplot(var$contrib, is.corr=FALSE)   
# Contributions of gene variables to PC1
fviz_contrib(res.pca, choice="var", axes = 1,top = 10, 
             fill = "lightgray", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45))
# Contributions of gene variables to PC2
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
quanti.correlation.dim1 = dim1$quanti.correlation 
names(quanti.correlation.dim1) = dim1$entrez
head(quanti.correlation.dim1)


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

##########################################################################
############ 12. Plot of normalized counts for a single gene  ############
##########################################################################
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
d## *************************************************************

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
}

