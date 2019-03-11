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

EIF.proteomics <- read.csv(file.path("project-data", 
                                     "proteomics.csv"), 
                           header = TRUE, 
                           sep = ",")
EIF.proteomics <- as.data.frame(EIF.proteomics)
EIF.proteomics <- EIF.proteomics[EIF.proteomics$gene_id != "", ]
EIF.proteomics <- droplevels(EIF.proteomics)
nlevels(EIF.proteomics$gene_id)
EIF.proteomics <- EIF.proteomics[!duplicated(EIF.proteomics$gene_id), ]
EIF.proteomics2 <- EIF.proteomics[ ,-1]
rownames(EIF.proteomics2) <- EIF.proteomics[ ,1]
guideData <- EIF.proteomics2
condition <- c("X60", "X40", "X20", "X10", "X5")
guideDesign <- data.frame(row.names = colnames(guideData),
                          condition = condition)
dds <- DESeqDataSetFromMatrix(countData = guideData,
                              colData   = guideDesign,
                              design    = ~ condition)
dds
head(assay(dds))

