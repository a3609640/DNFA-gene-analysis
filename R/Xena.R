## the following script perform PCA on RNA-Seq data of seven DNFA genes from SKCM amd GTEX dataset with R package "ggfortify".
library(cluster)
library(ggfortify)
library(ggplot2)
library(lfda)
library(readxl)

doAllXena <- function() {
DNFASKCMandGTEX <- readxl::read_excel(file.path("project-data", "7DNFASKCMandGTEX.xls"))

# eliminate unique row with "Solid Tissue Normal"
DNFASKCMandGTEX <- DNFASKCMandGTEX[-nrow(DNFASKCMandGTEX),]


View(DNFASKCMandGTEX)
str(DNFASKCMandGTEX)
DNFASKCMandGTEX <- as.data.frame(DNFASKCMandGTEX)
DNFASKCMandGTEX$sample_type <- as.factor(DNFASKCMandGTEX$sample_type)

DNFASKCMandGTEX$sample_type <- factor(DNFASKCMandGTEX$sample_type,
                                      levels = c("Normal Tissue",
                                                 "Primary Tumor",
                                                 "Metastatic"))
                                                 # "Solid Tissue Normal"))

levels(DNFASKCMandGTEX$sample_type)
length(DNFASKCMandGTEX)
list(DNFASKCMandGTEX)
# DNFASKCMandGTEX data table contain the RNA-Seq results for 7 DNFA genes,
# and 8th column shows the cancer category for each sample.
# data.frame for prcomp() function should only contains numeric values, so
# we remove the 8th column from DNFASKCMandGTEX data
df <- DNFASKCMandGTEX[-8]
plot1 <- ggplot2::autoplot(prcomp(df), data = DNFASKCMandGTEX, colour = 'sample_type')
print(plot1)
# we can color each data point in the plot according to their cancer category: colour = 'sample_type'
# Use loadings = TRUE, we draw eigenvector for each DNFA gene on the plot.
bp <- ggplot2::autoplot(prcomp(df), data = DNFASKCMandGTEX, colour = 'sample_type',
         loadings = TRUE, loadings.colour = 'black',
         loadings.label.vjust = -1,
         loadings.label = TRUE, loadings.label.size = 4) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = 'transparent', color = 'black', size = 1),
        axis.title = element_text(colour = "black", size = 12,
                                    face = "bold"),
        axis.text = element_text(colour = "black", size = 12,
                                 face = "bold"),
        legend.title = element_text(colour = "black", size = 12,
                                     face = "bold"),
        legend.text = element_text(colour = "black", size = 12,
                                 face = "bold", hjust = 1),
        legend.key = element_blank())

print(bp)

# use promp function only to draw the PCA plot for DNFA gene expression, similar plot as the past section.
PC <- prcomp(df)
PCi <- data.frame(PC$x, Species = DNFASKCMandGTEX$sample_type)
plot3 <- ggplot(PCi,aes(x = PC1, y = PC2, col = Species)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_color_manual(values = c("#FF1BB3","#A7FF5B","#99554D","#CC0000")) + #choose colors here
  theme_classic()
print(plot3)
# clara: Clustering Large Applications, we chose the number of clusters as "3".
# tumor tissues and healthy tissus samples are clearly seperated,
# However, metastatic and primary tumor samples can not be seperated by clustering method
# It is required that 0 < k < n where n is the number of factors. Here 4 categories of tissues
plot4 <- ggplot2::autoplot(clara(df, 2))
print(plot4)
# fanny: Fuzzy Analysis Clustering, to compute a fuzzy clustering of the data into k clusters.
plot5 <- ggplot2::autoplot(fanny(df, 2), frame = TRUE)
print(plot5)
# pam Partitioning Around Medoids, , a more robust version of K-means.
plot6 <- ggplot2::autoplot(pam(df, 2), frame = TRUE, frame.type = 'norm')
print(plot6)

# Dimensionality (r) is set to 7 because there are 7 DNFA genes to consider.

# Note, set knn = 1 and minObsPerLabel = 1 if there is a label in the data that only
# occurs once.  However, still expect errors; lfda has strong assumptions that there
# are multiple observations per label.
model <- lfda(DNFASKCMandGTEX[-8], DNFASKCMandGTEX[, 8], r = 7, metric = "plain",
              knn = 5)
plot7 <- ggplot2::autoplot(model, data = DNFASKCMandGTEX, frame = TRUE, frame.colour = 'sample_type')
print(plot7)

# A beta value of 0 indicates totally supervised learning; 1 is totally unsupervised.
model <- self(DNFASKCMandGTEX[-8], DNFASKCMandGTEX[, 8], beta = 0.0, r = 7, metric = "plain",
              minObsPerLabel = 5)
plot8 <- ggplot2::autoplot(model, data = DNFASKCMandGTEX, frame = TRUE, frame.colour = 'sample_type')
print(plot8)

model <- self(DNFASKCMandGTEX[-8], DNFASKCMandGTEX[, 8], beta = 1.0, r = 7, metric = "plain",
              minObsPerLabel = 5)
plot9 <- ggplot2::autoplot(model, data = DNFASKCMandGTEX, frame = TRUE, frame.colour = 'sample_type')
print(plot9)
}
doAllXena()
