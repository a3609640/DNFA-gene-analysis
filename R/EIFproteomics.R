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
library(gridExtra)
library(lattice)
library(Matrix)
library(org.Hs.eg.db)
library(pheatmap)
library(plyr)  # deliberately out of sorted order; must precede dplyr
library(RColorBrewer)
library(Rcpp)
library(RcppArmadillo)
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
colnames(guide.EIF.Data)<- EIF.condition
guide.EIF.Design <- data.frame(row.names = colnames(guide.EIF.Data),
                               condition = EIF.condition)


X4EGI_1_down <- read_csv("/home/wagner/suwu/Documents/translation/proteomics/4EGI-1 down.csv")
X4EGI_1_down <- data.frame(X4EGI_1_down)
X4EGI_1_down <- X4EGI_1_down
X4EGI_1_down <- X4EGI_1_down[!duplicated(X4EGI_1_down$GENE_SYMBOL), ]
rownames(X4EGI_1_down) <- X4EGI_1_down$GENE_SYMBOL
colnames(X4EGI_1_down)[colnames(X4EGI_1_down)=="GENE_ID"] <- "entrez"
EGI.down.genes <- subset(X4EGI_1_down, select=c("X4G.I_FoldChange"))
EIF.proteomics2$entrez = mapIds(org.Hs.eg.db,
                         keys      = row.names(EIF.proteomics2),
                         column    = "ENTREZID",
                         keytype   = "SYMBOL",
                         multiVals = "first")
EIF.proteomics2$name =   mapIds(org.Hs.eg.db,
                         keys      = row.names(EIF.proteomics2),
                         column    = "GENENAME",
                         keytype   = "SYMBOL",
                         multiVals = "first")


EIF.proteomics2$rn <- rownames(EIF.proteomics2)
EGI.down.genes$rn <- rownames(EGI.down.genes)
df <- join_all(list(EIF.proteomics2, X4EGI_1_down), by = 'entrez', type = 'full')


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


grid.arrange(p1, p2, p3, p4, nrow = 2)

df <- scale(guide.EIF.Data)
# Optimal number of clusters for k-means

my_data <- scale(guide.EIF.Data)
fviz_nbclust(my_data, kmeans, method = "gap_stat")

# 2. Compute k-means
km.res <- kmeans(scale(guide.EIF.Data), 1, nstart = 25)
# 3. Visualize
fviz_cluster(km.res, data = df,
             palette = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"),
             ggtheme = theme_minimal(),
             main = "Partitioning Clustering Plot"
)


#####################################
## 6. Principal component analysis ##
#####################################

res.pca <- PCA(t(guide.EIF.Data),  graph = FALSE)
# Extract eigenvalues/variances
get_eig(res.pca)
# Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 100))
# Extract the results for variables
var <- get_pca_var(res.pca)
var
# Coordinates of variables
head(var$coord)

# The function fviz_pca_ind() is used to visualize individuals
fviz_pca_ind(res.pca)
fviz_pca_ind(res.pca, col.ind = "cos2") +
  scale_color_gradient2(low      = "white", 
                        mid      = "blue", 
                        high     = "red", 
                        midpoint = 0.50)

.des_facto_get_title_font <- function() {
  return(element_text(face  = "bold", 
                      color = "black", 
                      size  = 18))
  }


# Graph of variables: default plot
# Make a biplot of individuals and variables :
pcaBiplot <- fviz_pca_biplot(res.pca,
                             geom = c("point", "text", "arrows"), 
                             select.var = list(contrib = 10),
                             alpha.var ="cos2") +
  geom_point(size = 3) +
  theme_bw()  +
  xlim(-40, 20) +
  ylim(-10, 15) +
  theme(text              = .des_facto_get_title_font(),
        axis.text         = .des_facto_get_title_font(),
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
        legend.position      = c(0,1),
        legend.justification = c(-0.05, 1.05)) +
  scale_color_brewer(palette = "Dark2")
print(pcaBiplot)
}




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
summary(dim1)
head(dim1,10)

res.desc$Dim.2
dim2 <- data.frame(res.desc$Dim.2)
summary(dim2)
head(dim2,10)

############################################################
## 3.2. Add Entrez IDs, gene symbols, and full gene names ##
############################################################
.addGeneIdentifiers <- function(EIF.Data) {
  EIF.Data <- data.frame(EIF.Data)
  columns(org.Hs.eg.db)
  EIF.Data$entrez = mapIds(org.Hs.eg.db,
                           keys      = row.names(EIF.Data),
                           column    = "ENTREZID",
                           keytype   = "SYMBOL",
                           multiVals = "first")
  EIF.Data$name =   mapIds(org.Hs.eg.db,
                           keys      = row.names(EIF.Data),
                           column    = "GENENAME",
                           keytype   = "SYMBOL",
                           multiVals = "first")
  summary(EIF.Data)
  head(EIF.Data, 10)
  return(EIF.Data)
}

dim1.entrez <- .addGeneIdentifiers(dim1)
dim2.entrez <- .addGeneIdentifiers(dim2)
######################################################
## 5. KEGG pathway analysis on dim1 from proteomics data ##
######################################################

.keggAnalysis <- function(dim1.entrez) {
  data(kegg.gs)
  lapply(kegg.gs[1:3],head)
## kegg.sets.hs is a named list of 229 elements
## Each element is a character vector of member gene Entrez IDs
## for a single KEGG pathway
  data(kegg.sets.hs)
## sigmet.idx.hs is an index of numbers of signaling and
## metabolic pathways in kegg.set.gs.
  data(sigmet.idx.hs)
## kegg.sets.hs[sigmet.idx.hs] gives you the “cleaner” gene sets of signaling
## and metabolic pathways only.
  kegg.sets.hs.sigmet <-  kegg.sets.hs[sigmet.idx.hs]
  head(kegg.sets.hs.sigmet, 3)
## Generate a named vector of fold changes, where the names of the values
## are the Entrez gene IDs, for the gage() function
  quanti.correlation <- dim1.entrez$quanti.correlation
  names(quanti.correlation) <- dim1.entrez$entrez
  head(quanti.correlation)
# Get KEGG pathway with both metabolism and signaling pathways
  kg.hsa <- kegg.gsets("hsa")
  kegg.sigmet.idx <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
  keggres.sigmet.idx <- gage(quanti.correlation,
                             gsets    = kegg.sigmet.idx,
                             same.dir = FALSE)
  lapply(keggres.sigmet.idx, head,20)
# Get KEGG pathway with only metabolism pathways
  kegg.met.idx <- kg.hsa$kg.sets[kg.hsa$met.idx]
  keggres.met.idx <- gage(quanti.correlation,
                          gsets    = kegg.met.idx,
                          same.dir = FALSE)
  lapply(keggres.met.idx, head,10)
# Get KEGG pathway with only signaling pathways
  kegg.sig.idx <- kg.hsa$kg.sets[kg.hsa$sig.idx]
  keggres.sig.idx <- gage(quanti.correlation,
                          gsets    = kegg.sig.idx,
                          same.dir = FALSE)
  lapply(keggres.sig.idx, head,10)
}

####################################################
## 6. GO pathway analysis on dim1 from proteomics data ##
####################################################
.pathwayAnalysis <- function(dim1.entrez) {
  # set up GO database
  data(go.gs)
  lapply(go.gs[1:3],head)
  go.hs <- go.gsets(species = "human")
  go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP] # “Biological Process”
  go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF] # “Molecular Function”
  go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC] # “Cellular Component”
  # GO pathway analysis with Biological Process
  quanti.correlation <- dim1.entrez$quanti.correlation
  names(quanti.correlation) <- dim1.entrez$entrez
  head(quanti.correlation)
  fc.go.bp.p <- gage(quanti.correlation,
                     gsets    = go.bp.gs,
                     same.dir = FALSE)
  lapply(fc.go.bp.p, head,10)
#  write.table(fc.go.bp.p$greater, file = "fc.go.bp.p.greater.txt", sep = "\t")
#  write.table(fc.go.bp.p$less, file = "fc.go.bp.p.less.txt", sep = "\t")
  
  # GO pathway analysis with Molecular Process
  fc.go.mf.p <- gage(quanti.correlation, 
                     gsets = go.mf.gs,
                     same.dir = FALSE)
  lapply(fc.go.mf.p, head,10)
  # GO pathway analysis with Cellular Process
  fc.go.cc.p <- gage(quanti.correlation, 
                     gsets = go.cc.gs,
                     same.dir = FALSE)
  lapply(fc.go.cc.p, head,10)
}

.keggAnalysis(dim1.entrez)

