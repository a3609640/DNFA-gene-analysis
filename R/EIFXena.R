## the following script perform PCA on RNA-Seq data of seven DNFA genes 
## from SKCM amd GTEX dataset with R package "ggfortify".
library(cluster)
library(factoextra)
library(ggfortify)
library(ggplot2)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(lfda)
library(readxl)
library(reshape2)


EIF.TCGA.GTEX <- readxl::read_excel(file.path("project-data",
                                                "EIFTCGAGTEX.xlsx"))
<<<<<<< HEAD
EIFTCGAGTEX$`_sample_type` <- as.factor(EIFTCGAGTEX$`_sample_type`)
EIFTCGAGTEX$`_study` <- as.factor(EIFTCGAGTEX$`_study`)
levels(EIFTCGAGTEX$`_sample_type`)
# 7DNFASKCMandGTEX.xls contains data downloaded from Xena
# TODO SW, will write a seperate script (with Xenatool) to download the data
# eliminate unique row with "Solid Tissue Normal"
EIFTCGAGTEX <- droplevels(EIFTCGAGTEX[!EIFTCGAGTEX$`_sample_type` == 'Solid Tissue Normal',])
EIFTCGAGTEX <- droplevels(EIFTCGAGTEX[!EIFTCGAGTEX$`_sample_type` == 'Cell Line',])
EIFGTEX <- EIFTCGAGTEX[EIFTCGAGTEX$`_study` == 'GTEX',]
EIFTCGA <- EIFTCGAGTEX[EIFTCGAGTEX$`_study` == 'TCGA',]
EIFGTEX <- EIFGTEX[, c(3,4,5,6)]
EIFTCGA <- EIFTCGA[, c(3,4,5,6)]
par(mfrow=c(1,2))
boxplot(EIFGTEX, main="EIF RNAseq in healthy tissues, n=8152")
boxplot(EIFTCGA, main="EIF RNAseq in tumor tissues, n=10531")
=======
EIF.TCGA.GTEX$`_sample_type` <- as.factor(EIF.TCGA.GTEX$`_sample_type`)
EIF.TCGA.GTEX$`_study` <- as.factor(EIF.TCGA.GTEX$`_study`)
levels(EIF.TCGA.GTEX$`_sample_type`)
EIF.TCGA.GTEX <- droplevels(EIF.TCGA.GTEX[!EIF.TCGA.GTEX$`_study` == 'TARGET',])
EIF.GTEX <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$`_study` == 'GTEX',]
EIF.TCGA.GTEX.long.form <- melt(EIF.TCGA.GTEX)
colnames(EIF.TCGA.GTEX.long.form) <- c("sample","study","sample.type",
                                     "primary.site","variable","value")
plotEIF(EIF.TCGA.GTEX.long.form)

EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$`_study` == 'TCGA',]
levels(EIF.TCGA$`_sample_type`)
tumor.type <- c("Metastatic","Primary Tumor",
                "Recurrent Tumor","Normal Tissue",
                "Solid Tissue Normal")
EIF.TCGA <- EIF.TCGA[EIF.TCGA$`_sample_type` %in% tumor.type,]
EIF.TCGA.long.form <- melt(EIF.TCGA)
colnames(EIF.TCGA.long.form) <- c("sample","study","sample.type",
                                  "primary.site","variable","value")
plotEIF(EIF.TCGA.long.form)

plotEIF <-  function (x) {
  name <- deparse(substitute(x))
  black_bold_tahoma_12 <- element_text(color  = "black", 
                                       face   = "bold",
                                       family = "Tahoma", 
                                       size   = 12)
  
  black_bold_tahoma_12_45 <- element_text(color  = "black",
                                          face   = "bold",
                                          family = "Tahoma", 
                                          size   = 12, 
                                          angle  = 45,
                                          hjust  = 1)
  ggplot(data = x,
         aes(x     = sample.type, 
             y     = value, 
             color = sample.type)) +
    facet_wrap( ~ variable) +
    geom_boxplot(alpha    = .01, 
                 size     = .75,
                 width    = .75,
                 position = position_dodge(width = .9)) +
    #    labs(title = paste0(name," n = 8555"),
    #         x     = "eIF4F complex components",
    #         y     = paste0("log2(value)")) +
    ylim(5, 25) +
    theme_bw() +
    theme(plot.title  = black_bold_tahoma_12,
          axis.title  = black_bold_tahoma_12,
          axis.text.x = black_bold_tahoma_12_45,
          axis.text.y = black_bold_tahoma_12,
          axis.line.x = element_line(color  = "black"),
          axis.line.y = element_line(color  = "black"),
          panel.grid  = element_blank(),
          strip.text  = black_bold_tahoma_12)
}


EIF.TCGA.long.form.EIF4A1 <- EIF.TCGA.long.form[EIF.TCGA.long.form$variable == "EIF4A1",]
EIF.TCGA.long.form.EIF4E <- EIF.TCGA.long.form[EIF.TCGA.long.form$variable == "EIF4E",]
EIF.TCGA.long.form.EIF4G1 <- EIF.TCGA.long.form[EIF.TCGA.long.form$variable == "EIF4G1",]
EIF.TCGA.long.form.EIF4EBP1 <- EIF.TCGA.long.form[EIF.TCGA.long.form$variable == "EIF4EBP1",]
my_comparison <- list(c("Metastatic", "Solid Tissue Normal"), 
                      c("Primary Tumor", "Solid Tissue Normal"), 
                      c("Recurrent Tumor", "Solid Tissue Normal"),
                      c("Metastatic", "Primary Tumor"),
                      c("Recurrent Tumor", "Primary Tumor"))


p1 <- plotEIF(EIF.TCGA.long.form.EIF4A1) +
  stat_compare_means(comparisons = my_comparison, method = "t.test")
p2 <- plotEIF(EIF.TCGA.long.form.EIF4E) +
  stat_compare_means(comparisons = my_comparison, method = "t.test")
p3 <- plotEIF(EIF.TCGA.long.form.EIF4G1) +
  stat_compare_means(comparisons = my_comparison, method = "t.test")
p4 <- plotEIF(EIF.TCGA.long.form.EIF4EBP1) +
  stat_compare_means(comparisons = my_comparison, method = "t.test")

grid.arrange(p1, p2, p3, p4, ncol = 2)
















>>>>>>> 9321f8f65304ea06f9845bdc4108446f64ece6b8

View(DNFASKCMandGTEX)
str(DNFASKCMandGTEX)
EIFTCGAGTEX <- as.data.frame(EIFTCGAGTEX)
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
plot1 <- ggplot2::autoplot(prcomp(df),
                           data   = DNFASKCMandGTEX,
                           colour = 'sample_type')
print(plot1)
.getTitleFont <- function() {
  return(element_text(face = "bold", color = "black", size = 12))
}
# we can color each data point in the plot according to their cancer category: 
# colour = 'sample_type'
# Use loadings = TRUE, we draw eigenvector for each DNFA gene on the plot.
bp <- ggplot2::autoplot(prcomp(df),
 data                   = DNFASKCMandGTEX,
 colour                 = 'sample_type', ## use colour, color gives error
 loadings               = TRUE,
 loadings.colour        = 'black',
 loadings.label.vjust   = -1,
 loadings.label         = TRUE,
 loadings.label.size    = 4) +
 theme(plot.background  = element_blank(),
       panel.background = element_rect(fill   = 'transparent',
                                       color  = 'black',
                                       size   = 1),
       axis.title       = .getTitleFont(),
       axis.text        = .getTitleFont(),
       legend.title     = .getTitleFont(),
       legend.text      = .getTitleFont(),
       legend.key       = element_blank())
print(bp)

# use promp function only to draw the PCA plot for DNFA gene expression, 
# similar plot as the past section.
PC <- prcomp(df)
PCi <- data.frame(PC$x, Species = DNFASKCMandGTEX$sample_type)
plot3 <- ggplot(PCi,aes(x = PC1, y = PC2, col = Species)) +
  geom_point(size = 3, alpha = 0.5) +
  # choose colors here
  scale_color_manual(values = c("#FF1BB3","#A7FF5B","#99554D","#CC0000")) + 
  theme_classic()
print(plot3)

## determine the optimal number of clusters K for k-means clustering:
fviz_nbclust(df, kmeans, method = "silhouette")
fviz_nbclust(df, pam, method = "silhouette")
# From the plot, the suggested optimal number of clusters is 2.
# clara: Clustering Large Applications is an extension to the PAM 
# (Partitioning Around Medoids). We chose the number of clusters as "2".
plot4 <- ggplot2::autoplot(clara(df, 2))
print(plot4)
# fanny: Fuzzy Analysis Clustering, 
# to compute a fuzzy clustering of the data into k clusters.
plot5 <- ggplot2::autoplot(fanny(df, 2), frame = TRUE)
print(plot5)
# pam: Partitioning Around Medoids, a more robust version of K-means.
plot6 <- ggplot2::autoplot(pam(df, 2), frame = TRUE, frame.type = 'norm')
print(plot6)
# Tumor tissues and healthy tissus samples are clearly seperated,
# Metastatic and primary tumor samples were not seperated by clustering method

# Local Fisher Discriminant Analysis (LFDA)
# Dimensionality (r) is set to 7 because there are 7 DNFA genes to consider.
# Note, set knn = 1 and minObsPerLabel = 1 if there is a label in the data 
# that only occurs once.  However, still expect errors; 
# lfda has strong assumptions that there are multiple observations per label.
model <- lfda(DNFASKCMandGTEX[-8],
              DNFASKCMandGTEX[, 8],
              r      = 7,
              metric = "plain",
              knn    = 5)
plot7 <- ggplot2::autoplot(model,
                           data         = DNFASKCMandGTEX,
                           frame        = TRUE,
                           frame.colour = 'sample_type')
print(plot7)

# A beta value of 0 indicates totally supervised learning; 
# A beta value of 1 indicates totally unsupervised.
model <- self(DNFASKCMandGTEX[-8],
              DNFASKCMandGTEX[, 8],
              beta           = 0.0,
              r              = 7,
              metric         = "plain",
              minObsPerLabel = 5)
plot8 <- ggplot2::autoplot(model,
                           data         = DNFASKCMandGTEX,
                           frame        = TRUE,
                           frame.colour = 'sample_type')
print(plot8)

model <- self(DNFASKCMandGTEX[-8],
              DNFASKCMandGTEX[, 8],
              beta           = 1.0,     # unsupervised learning
              r              = 7,
              metric         = "plain",
              minObsPerLabel = 5)
plot9 <- ggplot2::autoplot(model,
                           data         = DNFASKCMandGTEX,
                           frame        = TRUE,
                           frame.colour = 'sample_type')
print(plot9)
}
doAllXena()
