library(AnnotationDbi)
library(BiocStyle)
library(cluster)
library(colorspace)
library(ComplexHeatmap)
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
library(made4)
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
guide.EIF.Data <- EIF.proteomics2
EIF.condition <- c("60 uM", "40 uM", "20 uM", "10 uM", "5 uM")
guide.EIF.Design <- data.frame(row.names = colnames(guide.EIF.Data),
                               condition = EIF.condition)


dists <- dist(t(guide.EIF.Data))
plot(hclust(dists), labels = guide.EIF.Design$condition)
sampleDistMatrix <- as.matrix(dists)
rownames(sampleDistMatrix) <- paste("4EGi", EIF.condition, sep = "-")
colnames(sampleDistMatrix) <- paste("4EGi", EIF.condition, sep = "-")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = dists,
         clustering_distance_cols = dists,
         col                      = colors)


Heatmap(guide.EIF.Data)
EIFMatrix = data.matrix(guide.EIF.Data)
EIFMatrix = EIFMatrix[sample(nrow(EIFMatrix), nrow(EIFMatrix)), 
                      sample(ncol(EIFMatrix), ncol(EIFMatrix))]
Heatmap(EIFMatrix)



