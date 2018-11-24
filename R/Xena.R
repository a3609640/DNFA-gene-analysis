## the following script perform PCA on RNA-Seq data of seven DNFA genes from SKCM amd GTEX dataset with R package "ggfortify".  
library(readxl)
library(ggfortify)

DNFASKCMandGTEX <- read_excel("7DNFASKCMandGTEX.xls")
View(DNFASKCMandGTEX)
str(DNFASKCMandGTEX)
DNFASKCMandGTEX<-as.data.frame(DNFASKCMandGTEX)
DNFASKCMandGTEX$sample_type<-as.factor(DNFASKCMandGTEX$sample_type)
DNFASKCMandGTEX$sample_type <- factor(DNFASKCMandGTEX$sample_type, 
                                      levels = c("Normal Tissue", "Primary Tumor", "Metastatic", "Solid Tissue Normal"))

levels(DNFASKCMandGTEX$sample_type)
length(DNFASKCMandGTEX)
list(DNFASKCMandGTEX)
# DNFASKCMandGTEX data table contain the RNA-Seq results for 7 DNFA genes, and 8th columns shows the cancer category for each sample.
# data.frame for prcomp() function should only contains numeric values, so we remove the 8th column from DNFASKCMandGTEX data
df <- DNFASKCMandGTEX[c(1, 2, 3, 4, 5, 6, 7)]
autoplot(prcomp(df), data = DNFASKCMandGTEX, colour = 'sample_type')
# we can color each data point in the plot according to their cancer category: colour = 'sample_type'
# Use loadings = TRUE, we draw eigenvector for each DNFA gene on the plot.
bp = autoplot(prcomp(df), data = DNFASKCMandGTEX, colour = 'sample_type',
         loadings = TRUE, loadings.colour = 'black',
         loadings.label.vjust = -1,
         loadings.label = TRUE, loadings.label.size = 4)+
  theme(plot.background=element_blank(),
        panel.background=element_rect(fill='transparent',color='black',size=1),
        axis.title = element_text(colour="black", size=12,
                                    face="bold"),
        axis.text = element_text(colour="black", size=12,
                                 face="bold"),
        legend.title = element_text(colour="black", size=12,
                                     face="bold"),
        legend.text=element_text(colour="black", size=12,
                                 face="bold",hjust=1),
        legend.key=element_blank())

bp

# use promp function only to draw the PCA plot for DNFA gene expression, similar plot as the past section.
PC<-prcomp(df)
PCi<-data.frame(PC$x,Species=DNFASKCMandGTEX$sample_type)
ggplot(PCi,aes(x=PC1,y=PC2,col=Species))+
  geom_point(size=3,alpha=0.5)+ 
  scale_color_manual(values = c("#FF1BB3","#A7FF5B","#99554D","#CC0000"))+ #choose colors here
  theme_classic()


library(cluster)
# clara: Clustering Large Applications, we chose the number of clusters as "3". 
# tumor tissues and healthy tissus samples are clearly seperated,
# However, metastatic and primary tumor samples can not be seperated by clustering method
# It is required that 0 < k < n where n is the number of factors. Here 4 categories of tissues
autoplot(clara(DNFASKCMandGTEX[-7], 2))
# fanny: Fuzzy Analysis Clustering, to compute a fuzzy clustering of the data into k clusters.
autoplot(fanny(DNFASKCMandGTEX[-7], 2), frame = TRUE)
# pam Partitioning Around Medoids, , a more robust version of K-means.
autoplot(pam(DNFASKCMandGTEX[-7], 2), frame = TRUE, frame.type = 'norm')

library(lfda)
# Local Fisher Discriminant Analysis (LFDA)
model <- lfda(DNFASKCMandGTEX[,-7], DNFASKCMandTGEX[,7], r = 3, metric="plain")
autoplot(model, data = DNFASKCMandGTEX, frame = TRUE, frame.colour = 'sample_type')

model <- self(DNFASKCMandGTEX[-7], DNFASKCMandTGEX[, 7], beta = 0.1, r = 3, metric="plain")
autoplot(model, data = DNFASKCMandGTEX, frame = TRUE, frame.colour = 'sample_type')

model <- lfda(iris[-5], iris[, 5], r = 3, metric="plain")
autoplot(model, data = iris, frame = TRUE, frame.colour = 'Species')

model <- self(iris[-5], iris[, 5], beta = 0.1, r = 3, metric="plain")
autoplot(model, data = iris, frame = TRUE, frame.colour = 'Species')
