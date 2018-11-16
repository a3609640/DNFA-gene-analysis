library(plyr)
library(ggplot2) 
library(gplots)
library(grid)
library(data.table)

# singleRNAseq <- fread("GSE72056_melanoma_single_cell_revised_v2.txt")

get_lipogenesis_data <- function(data) {
  # -----------------------------------------------------------------------------------------------
  #select the genes of interest
  SREBF1 <- data[data$Cell=="SREBF1",] 
  SREBF2 <- data[data$Cell=="SREBF2",]
  FASN <- data[data$Cell=="FASN",] 
  SCD <- data[data$Cell=="SCD",] 
  ACACA <- data[data$Cell=="ACACA",]
  ACSS2 <- data[data$Cell=="ACSS2",]
  ACSL1 <- data[data$Cell=="ACSL1",]
  ACLY <- data[data$Cell=="ACLY",]
  HMGCR <- data[data$Cell=="HMGCR",]
  HMGCS1 <- data[data$Cell=="HMGCS1",]
  PPARGC1A <- data[data$Cell=="PPARGC1A",]
  PPARGC1B <- data[data$Cell=="PPARGC1B",]
  MITF <- data[data$Cell=="MITF",]
  AXL <- data[data$Cell=="AXL",]
  tumor <- data[data$Cell=="tumor",]
  malignant <- data[data$Cell=="malignant(1=no,2=yes,0=unresolved)",]
  
  # -----------------------------------------------------------------------------------------------
  #combine them into a new table "geneset"
  geneset <- rbind(tumor,malignant,SREBF1,SREBF2,FASN,SCD,ACACA,ACSS2,ACLY,ACSL1,HMGCR,HMGCS1,PPARGC1A,PPARGC1B,MITF,AXL)
  
  ## transpose the table
  # first remember the names
  n <- geneset$Cell
  # transpose all but the first column (name)
  transposed_geneset <- as.data.frame(t(geneset[,-1]))
  colnames(transposed_geneset) <- n
  renamed_geneset <- rename(
    transposed_geneset,
    c('malignant(1=no,2=yes,0=unresolved)'='malignancy'))
  renamed_geneset$malignancy <- factor(renamed_geneset$malignancy)
  renamed_geneset$malignancy <- revalue(renamed_geneset$malignancy,
                                        c("0" = "Unresolved",
                                          "1" = "Non-malignant",
                                          "2" = "Malignant"))
  return(renamed_geneset)
}

doAll4 <- function() {

lipogenesis_data <- get_lipogenesis_data(singleRNAseq)

# select RNA-seq data from malignant or non-malignant cells
totalgeneset <- subset(lipogenesis_data, malignancy != "Unresolved")
totalgeneset$tumor <- as.factor(totalgeneset$tumor)
levels(totalgeneset$tumor)
malignantgeneset <- subset(lipogenesis_data, malignancy == "Malignant")
nonmalignantgeneset <- subset(lipogenesis_data, malignancy == "Non-malignant")

#############################################################################################################################
#compare single cell SREBF1 expression level in malignant and non-malignant cells using boxplot graph within the same graph##
#############################################################################################################################
ggplot(totalgeneset, aes(x=factor(tumor), y=SREBF1)) + 
  geom_boxplot(size = 1, 
               outlier.colour=NA, 
               aes(fill = malignancy)) + 
  labs(x = "tumor samples", y = "SREBF1 mRNA counts") + 
  theme_bw()+ 
  theme(axis.title=element_text(face="bold",size=18,color="black"),
        axis.text=element_text(size=18,face="bold",color="black"),
        axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=18,face="bold",colour = 'black'),
        legend.position=c(0,1), 
        legend.justification=c(-0.1,1.1), 
        legend.key.height = unit(2.2, 'lines'),
        legend.background = element_rect(fill = "white")) + 
  guides(fill=guide_legend(title=NULL)) + 
  scale_y_continuous(breaks=seq(0, 10, 2))+
  scale_fill_discrete(labels=c("Nonmalignant cells", "Malignant cells"))

# -----------------------------------------------------------------------------------------------
ggplot(totalgeneset, aes(x=factor(tumor), y=FASN)) + 
  geom_boxplot(size = 1, 
               outlier.colour=NA, 
               aes(fill = malignancy)) + 
  labs(x = "tumor samples", y = "FASN mRNA counts") + 
  theme_bw()+ 
  theme(axis.title=element_text(face="bold",size=18,color="black"),
        axis.text=element_text(size=18,face="bold",color="black"),
        axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=18,face="bold",colour = 'black'),
        legend.position=c(0,1), 
        legend.justification=c(-0.1,1.1), 
        legend.key.height = unit(2.2, 'lines'),
        legend.background = element_rect(fill = "white")) + 
  guides(fill=guide_legend(title=NULL)) + 
  scale_y_continuous(breaks=seq(0, 10, 2))+
  scale_fill_discrete(labels=c("Nonmalignant cells", "Malignant cells"))

# -----------------------------------------------------------------------------------------------
ggplot(totalgeneset, aes(x=factor(tumor), y=SCD)) + 
  geom_boxplot(size = 1, 
               outlier.colour=NA, 
               aes(fill = malignancy)) + 
  labs(x = "tumor samples", y = "SCD mRNA counts") + 
  theme_bw()+ 
  theme(axis.title=element_text(face="bold",size=18,color="black"),
        axis.text=element_text(size=18,face="bold",color="black"),
        axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=18,face="bold",colour = 'black'),
        legend.position=c(0,1), 
        legend.justification=c(-0.1,1.1), 
        legend.key.height = unit(2.2, 'lines'),
        legend.background = element_rect(fill = "white")) + 
  guides(fill=guide_legend(title=NULL)) + 
  scale_y_continuous(breaks=seq(0, 10, 2))+
  scale_fill_discrete(labels=c("Nonmalignant cells", "Malignant cells"))

# -----------------------------------------------------------------------------------------------
ggplot(totalgeneset, aes(x=factor(tumor), y=ACACA)) + 
  geom_boxplot(size = 1, 
               outlier.colour=NA, 
               aes(fill = malignancy)) + 
  labs(x = "tumor samples", y = "ACACA mRNA counts") + 
  theme_bw()+ 
  theme(axis.title=element_text(face="bold",size=18,color="black"),
        axis.text=element_text(size=18,face="bold",color="black"),
        axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=18,face="bold",colour = 'black'),
        legend.position=c(0,1), 
        legend.justification=c(-0.1,1.1), 
        legend.key.height = unit(2.2, 'lines'),
        legend.background = element_rect(fill = "white")) + 
  guides(fill=guide_legend(title=NULL)) + 
  scale_y_continuous(breaks=seq(0, 10, 2))+
  scale_fill_discrete(labels=c("Nonmalignant cells", "Malignant cells"))

# -----------------------------------------------------------------------------------------------
ggplot(totalgeneset, aes(x=factor(tumor), y=SREBF2)) + 
  geom_boxplot(size = 1, 
               outlier.colour=NA, 
               aes(fill = malignancy)) + 
  labs(x = "tumor samples", y = "SREBF2 mRNA counts") + 
  theme_bw()+ 
  theme(axis.title=element_text(face="bold",size=18,color="black"),
        axis.text=element_text(size=18,face="bold",color="black"),
        axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=18,face="bold",colour = 'black'),
        legend.position=c(0,1), 
        legend.justification=c(-0.1,1.1), 
        legend.key.height = unit(2.2, 'lines'),
        legend.background = element_rect(fill = "white")) + 
  guides(fill=guide_legend(title=NULL)) + 
  scale_y_continuous(breaks=seq(0, 10, 2))+
  scale_fill_discrete(labels=c("Nonmalignant cells", "Malignant cells"))

# -----------------------------------------------------------------------------------------------
ggplot(totalgeneset, aes(x=factor(tumor), y=MITF)) + 
  geom_boxplot(size = 1, 
               outlier.colour=NA, 
               aes(fill = malignancy)) + 
  labs(x = "tumor samples", y = "MITF mRNA counts") + 
  theme_bw()+ 
  theme(axis.title=element_text(face="bold",size=18,color="black"),
        axis.text=element_text(size=18,face="bold",color="black"),
        axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=18,face="bold",colour = 'black'),
        legend.position=c(0,1), 
        legend.justification=c(-0.1,1.1), 
        legend.key.height = unit(2.2, 'lines'),
        legend.background = element_rect(fill = "white")) + 
  guides(fill=guide_legend(title=NULL)) + 
  scale_y_continuous(breaks=seq(0, 10, 2))+
  scale_fill_discrete(labels=c("Nonmalignant cells", "Malignant cells"))

# -----------------------------------------------------------------------------------------------
ggplot(totalgeneset, aes(x=factor(tumor), y=AXL)) + 
  geom_boxplot(size = 1, 
               outlier.colour=NA, 
               aes(fill = malignancy)) + 
  labs(x = "tumor samples", y = "AXL mRNA counts") + 
  theme_bw()+ 
  theme(axis.title=element_text(face="bold",size=18,color="black"),
        axis.text=element_text(size=18,face="bold",color="black"),
        axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=18,face="bold",colour = 'black'),
        legend.position=c(0,1), 
        legend.justification=c(-0.1,1.1), 
        legend.key.height = unit(2.2, 'lines'),
        legend.background = element_rect(fill = "white")) + 
  guides(fill=guide_legend(title=NULL)) + 
  scale_y_continuous(breaks=seq(0, 10, 2))+
  scale_fill_discrete(labels=c("Nonmalignant cells", "Malignant cells"))

}