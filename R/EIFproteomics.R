library(AnnotationDbi)
library(Biobase)
library(BiocGenerics)
library(BiocStyle)
library(cluster)    # clustering algorithms
library(colorspace)
library(ComplexHeatmap) # better heatmap generator
library(DESeq2)
library(dtwclust) # cluster time series with dynamic time warping
library(dplyr)
library(factoextra) # clustering algorithms & visualization
library(FactoMineR)
library(ggdendro) # dendrograms
library(ggplot2)
library(lattice)
library(Matrix)
library(org.Hs.eg.db)
library(pheatmap)
library(plyr)  # deliberately out of sorted order; must precede dplyr
library(Rcpp)
library(RcppArmadillo)
library(RColorBrewer)
library(S4Vectors)
library(SummarizedExperiment)
library(stats4)
library(stringr)
library(TCseq) # dose-dependent clustering
library(tidyverse)  # data manipulation
library(TSclust) # cluster time series
library(tseries) # bootstrap




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


# Dissimilarity matrix among protein targets
d <- dist(guide.EIF.Data, method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "ward" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)

# Dissimilarity matrix among 4EGI dosages
dists <- dist(t(guide.EIF.Data))
plot(hclust(dists), labels = guide.EIF.Design$condition)

# Compute with agnes
hc2 <- agnes(guide.EIF.Data, method = "complete")
# Agglomerative coefficient
hc2$ac
## [1] 0.8531583
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
# function to compute coefficient
ac <- function(x) {
  agnes(guide.EIF.Data, method = x)$ac
}
map_dbl(m, ac)

hc3 <- agnes(guide.EIF.Data, method = "ward")
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 





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


pc <- tsclust(guide.EIF.Data, type = "partitional", k = 20L, 
              distance = "dtw_basic", centroid = "pam", 
              seed = 3247L, trace = TRUE,
              args = tsclust_args(dist = list(window.size = 20L)))
plot(pc)
hc <- tsclust(guide.EIF.Data, type = "hierarchical", k = 20L, 
              distance = "sbd", trace = TRUE,
              control = hierarchical_control(method = "average"))
#> 
#> Calculating distance matrix...
#> Performing hierarchical clustering...
#> Extracting centroids...
#> 
#>  Elapsed time is 0.164 seconds.
plot(hc)


distance <- get_dist(guide.EIF.Data)
fviz_dist(distance, gradient = list(low = "#00AFBB", 
                                    mid = "white", 
                                    high = "#FC4E07"))
k2 <- kmeans(guide.EIF.Data, centers = 2, nstart = 25)
str(k2)
fviz_cluster(k2, data = guide.EIF.Data)

k3 <- kmeans(guide.EIF.Data, centers = 3, nstart = 25)
k4 <- kmeans(guide.EIF.Data, centers = 4, nstart = 25)
k5 <- kmeans(guide.EIF.Data, centers = 5, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", 
                   data = guide.EIF.Data) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point", 
                   data = guide.EIF.Data) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point", 
                   data = guide.EIF.Data) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point", 
                   data = guide.EIF.Data) + ggtitle("k = 5")

library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2)
