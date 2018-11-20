
library(readxl)
DNFASKCMandTGEX <- read_excel("7DNFASKCMandTGEX.xls")
View(DNFASKCMandTGEX)

str(DNFASKCMandTGEX)
DNFASKCMandTGEX<-as.data.frame(DNFASKCMandTGEX)
DNFASKCMandTGEX$sample_type<-as.factor(DNFASKCMandTGEX$sample_type)
DNFASKCMandTGEX$sample_type <- factor(DNFASKCMandTGEX$sample_type, 
                                      levels = c("Normal Tissue", "Primary Tumor", "Metastatic", "Solid Tissue Normal"))

levels(DNFASKCMandTGEX$sample_type)
length(DNFASKCMandTGEX)
list(DNFASKCMandTGEX)
  # conduct PCA on training dataset
library(ggfortify)
df <- DNFASKCMandTGEX[c(1, 2, 3, 4, 5, 6, 7)]
autoplot(prcomp(df))
# autoplot(prcomp(df), data = DNFASKCMandTGEX, colour = 'sample_type')

autoplot(prcomp(df), data = DNFASKCMandTGEX, colour = 'sample_type',
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

scale_fill_discrete(breaks=c("Normal Tissue", "Primary Tumor", "Metastatic", "Solid Tissue Normal"))+

PC<-prcomp(df)
PCi<-data.frame(PC$x,Species=DNFASKCMandTGEX$sample_type)
ggplot(PCi,aes(x=PC1,y=PC2,col=Species))+
  geom_point(size=3,alpha=0.5)+ #Size and alpha just for fun
  scale_color_manual(values = c("#FF1BB3","#A7FF5B","#99554D"))+ #your colors here
  theme_classic()







library(cluster)
autoplot(clara(DNFASKCMandTGEX[-7], 3))
autoplot(fanny(DNFASKCMandTGEX[-7], 5), frame = TRUE)
autoplot(pam(DNFASKCMandTGEX[-7], 3), frame = TRUE, frame.type = 'norm')

library(lfda)
# Local Fisher Discriminant Analysis (LFDA)
model <- lfda(DNFASKCMandTGEX[,-7], DNFASKCMandTGEX[,7], r = 3, metric="plain")
autoplot(model, data = DNFASKCMandTGEX, frame = TRUE, frame.colour = 'sample_type')

model <- self(DNFASKCMandTGEX[-7], DNFASKCMandTGEX[, 7], beta = 0.1, r = 3, metric="plain")
autoplot(model, data = DNFASKCMandTGEX, frame = TRUE, frame.colour = 'sample_type')

model <- lfda(iris[-5], iris[, 5], r = 3, metric="plain")
autoplot(model, data = iris, frame = TRUE, frame.colour = 'Species')

model <- self(iris[-5], iris[, 5], beta = 0.1, r = 3, metric="plain")
autoplot(model, data = iris, frame = TRUE, frame.colour = 'Species')
