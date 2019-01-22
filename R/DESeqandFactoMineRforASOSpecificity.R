## This R script will use the final processed RNA-Seq file from STAR alignment  
## (gene-counts-from-Makefile.csv) to perform differential gene expression analysis. 
## The following script uses DESeq and PCA analysis to compare 
## the specificity of SREBF1-ASOs to pooled siRNA of SREBF1 (positive control).

library(AnnotationDbi)
library(BiocStyle)
library(cluster)
library(colorspace)
library(DESeq2)
library(Biobase)
library(BiocGenerics)
library(plyr)  # deliberately out of sorted order; must precede dplyr
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
library(rmarkdown)
library(Rcpp)
library(RcppArmadillo)
library(RColorBrewer)
library(S4Vectors)
library(SummarizedExperiment)
library(stats4)
library(stringr)
library(survival)

.getTitleFont <- function() {
  return(element_text(face = "bold", color = "black", size = 18))
}

.readGeneCounts <- function() {
  # obtain the count table of the experiment directly 
  # from a pre-saved file: gene-counts.csv.
  # The RNA-seq was aligned to human reference genome Hg38 by STAR aligner
  # read processed RNA-seq read data from file testseq.csv.
  testseq <- read.csv(file.path("project-data", "gene-counts-from-Makefile.csv"))
  # Use the column one (Ensemble names) as columnn names.
  testseq <- data.frame(testseq[,-1], row.names = testseq[,1])
  # Remove the first four rows 
  # (N_unmapped,N_multimapping,N_noFeature and N_ambiguous)
  testseq <- data.frame(testseq[c(-1,-2,-3,-4),])
  # remove non-numeric 'symbol col' 25, leaving 4 col X 6 tests
  testseq <- testseq[-25]
  return(testseq)
}

##########################################################
## 2. Preparing count matrices from the RNA-seq results ##
##########################################################
.getGuideData <- function(testseq) {
  ## check the distribution of RNA-Seq reads
  par(mar = c(3,12,2,1))
  boxplot(testseq, outline = FALSE, horizontal = TRUE, las = 1)
  ## Remove rows in which there are no reads or nearly no reads
  guideData <- testseq[rowSums(testseq) > 1,]
  head(guideData)
  dim(guideData)
  ## check how the data distribution with boxplot 
  ## after removing rows with no read
  par(mar = c(3,12,2,1))
  boxplot(guideData, outline = FALSE, horizontal = TRUE, las = 1)
  return(guideData)
}

.getGuideDesign <- function(guideData, condition) {
  ## create a design for our "modelling"
  ## each sample contains four technical replicates
  return(data.frame(row.names = colnames(guideData),
                            condition = condition))
}

#####################################################
## 3. Construct DESeqDataSet from the count matrix ##
#####################################################
.getDDS <- function(guideData, guideDesign, condition) {
  ## Construct DESeqDataSet with the count matrix,
  ## countData, and the sample information, colData
  dds <- DESeqDataSetFromMatrix(countData = guideData,
                                colData   = guideDesign,
                                design    = ~ condition)
  dds
  head(assay(dds))
  return(dds)
}

##################################################
## 4. Standard differential expression analysis ##
##################################################
# DESeq function performs a default analysis through the steps:
# (1) estimation of size factor: estimateSizeFactors
# (2) estimation of dispersion: estimateDispersions
# (3) Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
.getDDSRES <- function(ddsDE) {
  ddsres <- results(ddsDE)
  # summary(ddsres)
  res <- data.frame(ddsres)
  return(ddsres)
}

##################################################
## 5. Count data transformations for clustering ##
##################################################
## The regularized log transform can be obtained using the rlog() function.
## Regularized log transform is to stabilize the variance of the data
## and to make its distribution roughly symmetric
## The default “blinds” the normalization to the design.
## This is very important so as to not bias the analyses (e.g. class discovery)
## The running times are shorter when using blind=FALSE
## and if the function DESeq has already been run,
## because then it is not necessary to re-estimate the dispersion values.
## The assay function is used to extract the matrix of normalized value
.makeHeatMap <- function(guideDesign, ddsDE) {
  vsd <- vst(ddsDE, blind = FALSE)
  rld <- rlog(ddsDE, blind = FALSE)
  # Hierarchical clustering using rlog transformation
  dists <- dist(t(assay(rld)))
  plot(hclust(dists), labels = guideDesign$condition)
  sampleDistMatrix <- as.matrix(dists)
  rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep = "-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = dists,
           clustering_distance_cols = dists,
           col                      = colors)
  return(rld)
}

#####################################
## 6. Principal component analysis ##
#####################################
.makePcaPlot <- function(rld) {
  ## number of top genes to use for principal components,
  ## selected by highest row variance, 500 by default
  pcaData <- plotPCA(rld, intgroup = c( "condition"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ## Print 2D PCA plot
  pcaPlot <- ggplot(
    data <- pcaData,
    aes_string(x     = "PC1",
               y     = "PC2",
               color = "condition")) +
    geom_point(size  = 3) +
    theme_bw() +
    xlim(-10, 6) +
    ylim(-6, 6) +
    theme(text                 = .getTitleFont(),
          axis.text            = .getTitleFont(),
          axis.line.x          = element_line(color = "black", size = 1),
          axis.line.y          = element_line(color = "black", size = 1),
          axis.ticks           = element_line(size  = 1),
          axis.ticks.length    = unit(.25, "cm"),
          panel.grid.major     = element_blank(),
          panel.grid.minor     = element_blank(),
          panel.border         = element_rect(color = "black", size = 1),
          panel.background     = element_blank(),
          legend.position      = c(0,0),
          legend.justification = c(-0.05,-0.05)) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance"))
  # Why 'print'?  See here:
  # https://stackoverflow.com/questions/26643852/ggplot-plots-in-scripts-do-not-display-in-rstudio
  print(pcaPlot)
  return(pcaData)
}

## Print 3D PCA plot
.plotPCA3D <- function(object,
                       intgroup   = "condition",
                       ntop       = 5000,
                       returnData = FALSE){
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
  d <- data.frame(PC1   = pca$x[, 1],
                  PC2   = pca$x[, 2],
                  PC3   = pca$x[, 3],
                  group = group,
                  intgroup.df,
                  name  = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:3]
    return(d)
  }
  message("Generating plotly plot")
  p <- plotly::plot_ly(data  = d,
                       x     = ~PC1,
                       y     = ~PC2,
                       z     = ~PC3,
                       color = group,
                       mode  = "markers",
                       type  = "scatter3d")
  print(p)
  return(p)
}

##########################################################
## 7. Add Entrez IDs, gene symbols, and full gene names ##
##########################################################
.addGeneIdentifiers <- function(res) {
  columns(org.Hs.eg.db)
  res$symbol = mapIds(org.Hs.eg.db,
                      keys      = row.names(res),
                      column    = "SYMBOL",
                      keytype   = "ENSEMBL",
                      multiVals = "first")
  res$entrez = mapIds(org.Hs.eg.db,
                      keys      = row.names(res),
                      column    = "ENTREZID",
                      keytype   = "ENSEMBL",
                      multiVals = "first")
  res$name =   mapIds(org.Hs.eg.db,
                      keys      = row.names(res),
                      column    = "GENENAME",
                      keytype   = "ENSEMBL",
                      multiVals = "first")
  summary(res)
  head(res, 10)
  return(res)
}

######################################################
## 8. Annotate genes enriched in each PCA component ##
######################################################
## scale.unit : a logical value. If TRUE, the data are scaled to unit variance 
## before the analysis.
## This standardization to the same scale avoids some variables to 
## become dominant just because of their large measurement units.
## We used FAlSE for scale.unit because rld has been run with DESEQ function before.
## ncp : number of dimensions kept in the final results.
## graph : a logical value. If TRUE a graph is displayed.
.annotateRld <- function(rld, guideDesign) {
  head(assay(rld))
  assayrld <- assay(rld)
  Pvars <- rowVars(assayrld)
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(500,
                                                        length(Pvars)))]
  columns(org.Hs.eg.db)
  row.names(assayrld) <- mapIds(org.Hs.eg.db,
                                keys      = row.names(assayrld),
                                column    = "SYMBOL",
                                keytype   = "ENSEMBL",
                                multiVals = "first")
  assayrld <- data.frame(t(assayrld[select, ]))
  assayrld$condition = guideDesign$condition
  con <- c("Mock-1", "Mock-2", "Mock-3", "Mock-4",
          "siNegative-1", "siNegative-2", "siNegative-3", "siNegative-4",
          "siSREBF1-1", "siSREBF1-2", "siSREBF1-3", "siSREBF1-4",
          "ASO-Neg-1", "ASO-Neg-2", "ASO-Neg-3", "ASO-Neg-4",
          "ASO-1-1", "ASO-1-2","ASO-1-3","ASO-1-4",
          "ASO-4-1", "ASO-4-2", "ASO-4-3", "ASO-4-4")
  rownames(assayrld) = con
  return(assayrld)
}
# The variable Species (index = 501) is removed
# before PCA analysis
## # Compute PCA with ncp = 3, to keep only the first three principal components
.makeAnnotatedPcaPlot <- function(assayrld) {
  return(PCA(assayrld[,-501], scale.unit = FALSE, ncp = 2,graph = TRUE))
}

######################################################
## 9. Extract variances in each principal component ##
######################################################
## Eigenvalues correspond to the amount of the variation explained 
## by each principal component (PC).
## Eigenvalues are large for the first PC and small for the subsequent PCs.
.makeBiplot <- function(assayrld, res.pca, pcaData) {
  eigenvalues <- res.pca$eig
  head(eigenvalues[, 1:2])
  eigen <- eigenvalues[1:10,]
  # Make a scree plot using base graphics:
  # Scree plot is a graph of eigenvalues/variances associated with components.
  barplot(eigen[, 2],
          names.arg = 1:nrow(eigen),
          main      = "Variances",
          xlab      = "Principal Components",
          ylab      = "Percentage of variances",
          col       = "steelblue")
  lines(x = 1:nrow(eigen), eigen[, 2],
        type = "b", pch = 19, col = "red")

  # plot biplot graph with the top six contributing genes to PCA from RNA-Seq
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaBiplot <-
    fviz_pca_biplot(res.pca,
                    select.var = list(contrib = 6),
                    # select.var = list(contrib = 0.6),
                    col.var    = "red",
                    label      = "var",
                    habillage  = assayrld$condition) +
      geom_point(size = 3,
         aes(color = factor(assayrld$condition))) +
          theme_bw() +
          xlim(-8, 4) +
          ylim(-5, 5) +
          theme(text              = .getTitleFont(),
                axis.text         = .getTitleFont(),
                axis.line.x       = element_line(color = "black",
                                                 size  = 1),
                axis.line.y       = element_line(color = "black",
                                                 size  = 1),
                axis.ticks        = element_line(size  = 1),
                axis.ticks.length = unit(.25, "cm"),
                panel.grid.major  = element_blank(),
                panel.grid.minor  = element_blank(),
                panel.border      = element_rect(color = "black",
                                                 size  = 1),
                panel.background  = element_blank(),
                legend.text       = element_text(color = "black",
                                                 size  = 18,
                                                 face  = "bold"),
            legend.position = c(0,1),
            legend.justification = c(-0.05, 1.05)) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance"))
  print(pcaBiplot)
}

#########################################################
## 10. Hierarchical Clustering on Principal Components ##
#########################################################
## Compute hierarchical clustering: Hierarchical clustering is performed 
## using the Ward’s criterion on the selected principal components.
## Ward criterion is used in the hierarchical clustering 
## because it is based on the multidimensional variance 
## like principal component analysis.
.makeHierarchicalCluster <- function(res.pca, pcaData) {
  res.hcpc <- HCPC(res.pca, graph = FALSE)
  plot1 <-
    fviz_dend(res.hcpc,
              cex                 = 0.7,     # Label size
              palette             = "jco",   # Color palette see ?ggpubr::ggpar
              rect                = TRUE,
              rect_fill           = TRUE,    # Add rectangle around groups
              rect_border         = "jco",   # Rectangle color
              labels_track_height = 0.8      # Augment the room for labels
    )
  plot2 <-
    fviz_cluster(res.hcpc,
                 repel            = TRUE,    # Avoid label overlapping
                 show.clust.cent  = TRUE,    # Show cluster centers
                 palette          = "jco",   # Color palette see ?ggpubr::ggpar
                 ggtheme          = theme_minimal(),
                 main             = "Factor map"
    )
  print(plot1)
  print(plot2)
  # Hierarchical Clustering
  # compute dissimilarity matrix for all data
  eu.d <- dist(pcaData, method  = "euclidean")
  # Hierarchical clustering using Ward's method
  res.hc <- hclust(eu.d, method = "ward.D2" )
  # Cut tree into 4 groups
  grp <- cutree(res.hc, k = 2)
  # Visualize
  plot(res.hc, cex = 0.6) # plot tree
  rect.hclust(res.hc, k = 2, border = c("yellow","blue")) # add rectangle
}

########################################################
## 11. Top contributing gene variables to PC1 and PC2 ##
########################################################
.findTopPrincipalComponentContributors <- function(res.pca) {
  ## Contributions of variables to PCs
  head(res.pca$var$contrib,10)
  head(res.pca$var$cos2, 10)
  var <- get_pca_var(res.pca)
  var
  ## library("corrplot")
  ## corrplot(var$contrib, is.corr=FALSE)
  ## Contributions of gene variables to PC1
  pc1Plot <- fviz_contrib(res.pca,
                          choice = "var",
                          axes   = 1,
                          top    = 10,
                          fill   = "lightgray",
                          color  = "black") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45))
  ## Contributions of gene variables to PC2
  pc2Plot <- fviz_contrib(res.pca,
                          choice = "var",
                          axes   = 2,
                          top    = 10,
                          fill   = "lightgray",
                          color  = "black") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45))
  print(pc1Plot)
  print(pc2Plot)
  ## Identify the most correlated variables with a given principal component
  res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)
  ## Description of dimension 1
  head(res.desc, 10)
  res.desc$Dim.1
  dim1 <- data.frame(res.desc$Dim.1)
  columns(org.Hs.eg.db)
  ## TODO(dlroxe): Figure out why it is necessary to comment out
  ## gene identifier assignment, here and for Dim.2 below.
  ## head(dim1,10)
  ## dim1 <- .addGeneIdentifiers(dim1)
  summary(dim1)
  head(dim1,10)
  quanti.correlation.dim1        = dim1$quanti.correlation
  names(quanti.correlation.dim1) = dim1$entrez
  head(quanti.correlation.dim1)
  res.desc$Dim.2
  dim2 <- data.frame(res.desc$Dim.2)
  columns(org.Hs.eg.db)
  ## dim1 <- .addGeneIdentifiers(dim1)
  summary(dim2)
  head(dim2,10)
}

######################################################
## 12. Plot of normalized counts for a single gene  ##
######################################################
# plotcount: "normalized" whether the counts should be normalized
# by size factor (default is TRUE)
# plotcount: "transform" whether to present log2 counts (TRUE)
# or to present the counts on the log scale (FALSE, default)
# re-arrange x-ase according to the following order:
# "Mock","siNeg","siBF1","ASO-neg","ASO-1","ASO-4"
# theme_bw() removes background color in the graph,
# guides(fill=FALSE) removes legends

.getCountsAndConditions <- function(dds, ensgID) {
  geneCounts <- plotCounts(
    dds, gene  = ensgID,
    intgroup   = "condition",
    normalized = TRUE,
    transform  = TRUE,
    returnData = TRUE)
  geneCounts$condition <- factor(geneCounts$condition,
    levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))
  return(geneCounts)
}

.makeGeneCountPlot <- function(countsAndConditions, title, ylim1, ylim2) {
  normalizedGeneCountTheme <-
    theme_bw() +
    theme(text              = .getTitleFont(),
          axis.text         = .getTitleFont(),
          axis.text.x       = element_text(angle = 45,      hjust = 1),
          axis.line.x       = element_line(color = "black", size  = 1),
          axis.line.y       = element_line(color = "black", size  = 1),
          axis.ticks        = element_line(size  = 1),
          axis.ticks.length = unit(.25, "cm"),
          panel.grid.major  = element_blank(),
          panel.grid.minor  = element_blank(),
          panel.border      = element_rect(color = "black", size  = 1),
          panel.background  = element_blank())
  plot <-
    ggplot(
      countsAndConditions,
      aes(x = condition, y = log2(count), fill = condition)) +
    geom_boxplot() +
    ylim(ylim1, ylim2) +
    guides(fill = FALSE) +
    normalizedGeneCountTheme +
    labs(title = title, x = " ", y = "log2(read counts)")
  print(plot)
}

analyze_aso_specificity <- function() {
  testseq <- .readGeneCounts()
  guideData <- .getGuideData(testseq)
  condition <- c(rep("Mock",4),rep("siNegative",4),rep("siSREBF1",4),
                 rep("ASO-Neg",4),rep("ASO-1",4),rep("ASO-4",4))
  guideDesign <- .getGuideDesign(guideData, condition)
  dds <- .getDDS(guideData, guideDesign, condition)
  ddsDE <- DESeq(dds)
  res <- .getDDSRES(ddsDE)
  rld <- .makeHeatMap(guideDesign, ddsDE)

  pcaData <- .makePcaPlot(rld)
  .plotPCA3D(rld, intgroup = "condition", ntop = 5000, returnData = FALSE)

  res <- .addGeneIdentifiers(res)
  assayrld <- .annotateRld(rld, guideDesign)
  res.pca <- .makeAnnotatedPcaPlot(assayrld)

  .makeBiplot(assayrld, res.pca, pcaData)

  .makeHierarchicalCluster(res.pca, pcaData)

  .findTopPrincipalComponentContributors(res.pca)


  # TODO(dlroxe): Use a function to derive "ENSG.." from
  # gene name.  Then the following repeated calls can be
  # collapsed into a single function that iterates over
  # a list of SREBF1, SCD, etc.
  SREBF1 <- .getCountsAndConditions(dds, "ENSG00000072310")
  .makeGeneCountPlot(SREBF1, "SREBF1", 8, 10.5)

  SCD <- .getCountsAndConditions(dds, "ENSG00000099194")
  .makeGeneCountPlot(SCD, "SCD", 13, 14.5)

  FASN <- .getCountsAndConditions(dds, "ENSG00000169710")
  .makeGeneCountPlot(FASN, "FASN", 12, 13.25)

  ACACA <- .getCountsAndConditions(dds, "ENSG00000278540")
  .makeGeneCountPlot(ACACA, "ACACA", 9.7, 10.5)

  ACSL1 <- .getCountsAndConditions(dds, "ENSG00000169710")
  .makeGeneCountPlot(ACSL1, "ACSL1", 11.5, 13.5)

  BRAF <- .getCountsAndConditions(dds, "ENSG00000157764")
  .makeGeneCountPlot(BRAF, "BRAF", 6, 9)

  SPIRE1 <- .getCountsAndConditions(dds, "ENSG00000134278")
  .makeGeneCountPlot(SPIRE1, "SPIRE1", 4, 9)

  USP9X <- .getCountsAndConditions(dds, "ENSG00000124486")
  .makeGeneCountPlot(USP9X, "USP9X", 7, 10)

  INSIG1 <- .getCountsAndConditions(dds, "ENSG00000186480")
  .makeGeneCountPlot(INSIG1, "INSIG1", 10, 12)

  YY1 <- .getCountsAndConditions(dds, "ENSG00000100811")
  .makeGeneCountPlot(YY1, "YY1", 8, 10.5)

  FABP7 <- .getCountsAndConditions(dds, "ENSG00000164434")
  .makeGeneCountPlot(FABP7, "FABP7", 6, 9)

  MITF <- .getCountsAndConditions(dds, "ENSG00000187098")
  .makeGeneCountPlot(MITF, "MITF", 10, 12)

  SREBF2 <- .getCountsAndConditions(dds, "ENSG00000198911")
  .makeGeneCountPlot(SREBF2, "SREBF2", 10, 12)

  HMGCR <- .getCountsAndConditions(dds, "ENSG00000113161")
  .makeGeneCountPlot(HMGCR, "HMGCR", 10, 12)

  NRAS <- .getCountsAndConditions(dds, "ENSG00000213281")
  .makeGeneCountPlot(NRAS, "NRAS", 7, 9)
}

# if (interactive()) {
#  analyze_aso_specificity()
# }

analyze_aso_specificity()
