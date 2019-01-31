library(data.table)
library(stringr)
library(ggplot2)
setwd("~/CCLE")

CCLE.expression <- fread("CCLE_Expression_Entrez_2012-09-29.csv",showProgress = T)
CCLE.anotation <- fread("CCLE_sample_info_file_2012-10-18.csv",showProgress = T)
# CCLE.anotation <- CCLE.anotation[,c(1,5)]
CCLE.drug <- fread("CCLE_NP24.2009_Drug_data_2015.02.24.csv",showProgress = T)
CCLE.RNAseq <- fread("CCLE_RNAseq_rsem_genes_tpm_20180929.txt",showProgress = T)
CCLE.gene.counts <- fread("CCLE_RNAseq_genes_counts_20180929.gct",showProgress = T)



############################################################################################################################################################################################################
############################################################################################################################################################################################################
## combine SREBP1 gene expression with drug sensitivity data

plotGene <- function(gene.name) { 
      EIF.gene <- c("EIF4A1","EIF4B","EIF4E","EIF4G1","EIF4EBP1","RPS6KB1","MYC")
      gene.object <- CCLE.gene.counts[CCLE.gene.counts$Description %in% EIF.gene,]
      gene.object <- gene.object[, -1]
      gene.object <- t(gene.object)
      class(gene.object)
      str(gene.object)
      ## SREBP1 is generated as a matrix somehow, and it needs to be converted into data frame
      gene.object <- as.data.frame(gene.object)
      setDT(gene.object, keep.rownames = TRUE)[]
      gene.object <- merge(gene.object, CCLE.anotation,by = 'CCLE name')
      class(gene.object$"Site Primary")
      class(gene.object$"gene.name")
      gene.object <- gene.object[ ,c(1,2,3,6)]
      gene.object$gene.name <- as.numeric(levels(gene.object$gene.name))[gene.object$gene.name]
      names(gene.object) <- c("CCLE.name", "gene.name", 
                              "cell.line.name","Primary.site")
      mean <- within(gene.object, 
                     Primary.site <- reorder(Primary.site, 
                                             log2(gene.name), median))
      ggplot(mean, aes(x = Primary.site, 
                       y = log2(gene.name))) + 
        geom_boxplot() + 
        theme_bw() + 
        labs(x = "Primary Site (CCLE)",
             y = paste("CCLE RNAseq gene expression of ", gene.name)) +
        theme(axis.title   = element_text(face  = "bold",
                                          size  = 12,
                                          color = "black"),
              axis.text    = element_text(size  = 12,
                                          angle = 90, 
                                          hjust = 1, 
                                          face  = "bold",
                                          color = "black"),
              axis.line.x  = element_line(color = "black"),
              axis.line.y  = element_line(color = "black"),
              panel.grid   = element_blank(),
              strip.text   = element_text(face  = "bold", 
                                        size    = 12, 
                                        color   = "black"),
              legend.text  = element_text(size  = 12,
                                         face   = "bold",
                                         colour = 'black'),
              legend.title = element_text(face  = "bold", 
                                          size  = 12, 
                                          color = "black"),
              legend.justification = c(1,1)) + 
        guides(fill = guide_legend(title = NULL))
}


plotGene(gene.name = "MITF")
plotGene(gene.name = "EIF4E")
plotGene(gene.name = "EIF4G1")
plotGene(gene.name = "EIF4EBP1")
plotGene(gene.name = "RPS6KB1")
plotGene(gene.name = "MYC")

SREBF1expression$SREBF1 <- as.numeric(as.character(SREBF1expression$SREBF1))
## somehow the numbers of SREBF1 columns are all changed into character 
SREBF1skin <- subset(SREBF1expression, `Site Primary`=="skin")
write.csv(SREBF1skin,"SREBF1skin.csv")

############################################################################################################################################################################################################
############################################################################################################################################################################################################
## save cell line names
n <- CCLE.expression$Description
## Transposing a dataframe CCLE.expression but removing the first column (microarry index for genes)
CCLE.expression.t <- as.data.frame(t(CCLE.expression[,-1:-2]))
##  use the cell line names to rename the colunms of the new dataset
colnames(CCLE.expression.t) <- n
## select only the skin cancer cell lines from the expression dataset
skin.cell.lines <- subset(CCLE.anotation, `Site Primary`=="skin")
CCLE.expression.skin <- CCLE.expression.t[skin.cell.lines$'CCLE name', ]
CCLE.expression.skin.lipo <- CCLE.expression.skin[ ,c("SCD", "FASN", "SREBF1","HMGCS1","HMGCR","SREBF2","MITF")]
colnames(ccle.genomics.lipo) <- c("SCD", "FASN", "SREBF1","HMGCS1","HMGCR","SREBF2","MITF")

       


plotDrug <- function(goi) { 
  expression_drugdata_skin <- subset(expression_drugdata, `Tissue.descriptor.1`=="skin")
  expression_drugdata_skin_goi <- subset(expression_drugdata_skin, `TARGET_PATHWAY`==goi)
  ggplot(expression_drugdata_skin_goi, aes(x=log2(SCD), y=log2(RMSE), colour=factor(TARGET))) + geom_point()+ facet_wrap(~TARGET)+ theme_bw()+ 
    labs(x = "SCD expression log2(RMA)", y = "log2(RMSE)") +   
    theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 12, colour = "black"),
          legend.text = element_text(size=12,face="bold",colour = 'black'),
          legend.title = element_text(face = "bold", size = 12, colour = "black"),
          legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))
}

#plotDrug(goi="ERK MAPK signaling")
## use loop and funtion!!
x <- levels(expression_drugdata$TARGET_PATHWAY)
lapply(x, plotDrug)

##m <- CorDrugResult$estimate > 0.25


############################################################################################################################################################################################################


drug_data_short <-drug_data[,c(1,2,3,4,13)]
colnames(drug_data_short)[1] <- "CCLE name"
SREBF1expression_drug_data <- merge(SREBF1expression,drug_data_short,by='CCLE name')
class(SREBF1expression_drug_data$Compound)
as.factor(SREBF1expression_drug_data$Compound)
class(SREBF1expression_drug_data$SREBF1)
class(SREBF1expression_drug_data$ActArea)
SREBF1expression_drug_data$SREBF1 <- as.numeric(as.character(SREBF1expression_drug_data$SREBF1))
SREBF1expression_drug_data$ActArea <- as.numeric(as.character(SREBF1expression_drug_data$ActArea))

ggplot(SREBF1expression_drug_data, aes(x=SREBF1, y=ActArea, colour=factor(Compound))) + geom_point()+ facet_wrap(~Compound)+ theme_bw()+ 
  labs(x = "SREBF1 expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))

class(SREBF1expression_drug_data$Compound)
SREBF1expression_drug_data$Compound <- as.factor(SREBF1expression_drug_data$Compound)
SREBF1expression_drug_data_Paclitaxel <- subset(SREBF1expression_drug_data, `Compound`=="Paclitaxel")
cor.test(SREBF1expression_drug_data_Paclitaxel$SREBF1, SREBF1expression_drug_data_Paclitaxel$ActArea, method = "pearson")

SREBF1expression_drug_data_17AAG <- subset(SREBF1expression_drug_data, `Compound`=="17-AAG")
cor.test(SREBF1expression_drug_data_17AAG$SREBF1, SREBF1expression_drug_data_17AAG$ActArea, method = "pearson")


ggplot(SREBF1expression_drug_data, aes(x=SREBF1, y=ActArea, colour=factor(Target))) + geom_point()+ facet_wrap(~Target)+ theme_bw()+ 
  labs(x = "SREBF1 expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))


SREBF1expression_drug_data_skin <- subset(SREBF1expression_drug_data, `Site Primary`=="skin")
ggplot(SREBF1expression_drug_data_skin, aes(x=SREBF1, y=ActArea, colour=factor(Compound))) + geom_point()+ facet_wrap(~Compound)+ theme_bw()+ 
  labs(x = "SREBF1 expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))


class(SREBF1expression_drug_data_skin$Compound)
SREBF1expression_drug_data_skin$Compound <- as.factor(SREBF1expression_drug_data_skin$Compound)
SREBF1expression_drug_data_skin_MEK <- subset(SREBF1expression_drug_data_skin, `Compound`=="AZD6244")
cor.test(SREBF1expression_drug_data_skin_MEK$SREBF1, SREBF1expression_drug_data_skin_MEK$ActArea, method = "pearson")

SREBF1expression_drug_data_skin_MEK2 <- subset(SREBF1expression_drug_data_skin, `Compound`=="PD-0325901")
cor.test(SREBF1expression_drug_data_skin_MEK2$SREBF1, SREBF1expression_drug_data_skin_MEK2$ActArea, method = "pearson")

SREBF1expression_drug_data_skin_MEK2 <- subset(SREBF1expression_drug_data_skin, `Target`=="MEK")
cor.test(SREBF1expression_drug_data_skin_MEK2$SREBF1, SREBF1expression_drug_data_skin_MEK2$ActArea, method = "pearson")

SREBF1expression_drug_data_skin_HSP90 <- subset(SREBF1expression_drug_data_skin, `Compound`=="17-AAG")
cor.test(SREBF1expression_drug_data_skin_HSP90$SREBF1, SREBF1expression_drug_data_skin_HSP90$ActArea, method = "pearson")

SREBF1expression_drug_data_skin_Paclitaxel <- subset(SREBF1expression_drug_data_skin, `Compound`=="Paclitaxel")
cor.test(SREBF1expression_drug_data_skin_Paclitaxel$SREBF1, SREBF1expression_drug_data_skin_Paclitaxel$ActArea, method = "pearson")

ggplot(SREBF1expression_drug_data_skin, aes(x=SREBF1, y=ActArea, colour=factor(Target))) + geom_point()+ facet_wrap(~Target)+ theme_bw()+ 
  labs(x = "SREBF1 expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))

############################################################################################################################################################################################################
############################################################################################################################################################################################################
SREBF2 <- geneexpression[geneexpression$Description=="SREBF2",]
SREBF2<-t(SREBF2)
str(SREBF2)
## SREBP1 is generated as a matrix somehow, and it needs to be converted into data frame
SREBF2 <- as.data.frame(SREBF2)
setDT(SREBF2, keep.rownames = TRUE)[]
colnames(SREBF2) <- c("CCLE name","SREBF2")
#A one line option is: df$names<-rownames(df)
SREBF2expression <- merge(SREBF2,Annotations,by='CCLE name')
## somehow the numbers of SCD columns are all changed into character 
SREBF2skin <- subset(SREBF2expression, `Site Primary`=="skin")
write.csv(SREBF2skin,"SREBF2skin.csv")

SREBF2expression_drug_data <- merge(SREBF2expression,drug_data_short,by='CCLE name')
class(SREBF2expression_drug_data$Compound)
as.factor(SREBF2expression_drug_data$Compound)
class(SREBF2expression_drug_data$SREBF2)
class(SREBF2expression_drug_data$ActArea)
SREBF2expression_drug_data$SREBF2 <- as.numeric(as.character(SREBF2expression_drug_data$SREBF2))
ggplot(SREBF2expression_drug_data, aes(x=SREBF2, y=ActArea, colour=factor(Compound))) + geom_point()+ facet_wrap(~Compound)+ theme_bw()+ 
  labs(x = "SREBF2 expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))

ggplot(SREBF2expression_drug_data, aes(x=SREBF2, y=ActArea, colour=factor(Target))) + geom_point()+ facet_wrap(~Target)+ theme_bw()+ 
  labs(x = "SREBF2 expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))


SREBF2expression_drug_data_skin <- subset(SREBF2expression_drug_data, `Site Primary`=="skin")
ggplot(SREBF2expression_drug_data_skin, aes(x=SREBF2, y=ActArea, colour=factor(Compound))) + geom_point()+ facet_wrap(~Compound)+ theme_bw()+ 
  labs(x = "SREBF2 expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))


class(SREBF2expression_drug_data_skin$Compound)
SREBF2expression_drug_data_skin$Compound <- as.factor(SREBF2expression_drug_data_skin$Compound)
SREBF2expression_drug_data_skin_MEK <- subset(SREBF2expression_drug_data_skin, `Compound`=="AZD6244")
cor.test(SREBF2expression_drug_data_skin_MEK$SREBF2, SREBF2expression_drug_data_skin_MEK$ActArea, method = "pearson")

SREBF2expression_drug_data_skin_PD0325901 <- subset(SREBF2expression_drug_data_skin, `Compound`=="PD-0325901")
cor.test(SREBF2expression_drug_data_skin_PD0325901$SREBF2, SREBF2expression_drug_data_skin_PD0325901$ActArea, method = "pearson")

SREBF2expression_drug_data_skin_MEK2 <- subset(SREBF2expression_drug_data_skin, `Target`=="MEK")
cor.test(SREBF2expression_drug_data_skin_MEK2$SREBF2, SREBF2expression_drug_data_skin_MEK2$ActArea, method = "pearson")

SREBF2expression_drug_data_skin_HSP90 <- subset(SREBF2expression_drug_data_skin, `Compound`=="17-AAG")
cor.test(SREBF2expression_drug_data_skin_HSP90$SREBF2, SREBF2expression_drug_data_skin_HSP90$ActArea, method = "pearson")



ggplot(SREBF2expression_drug_data_skin, aes(x=SREBF2, y=ActArea, colour=factor(Target))) + geom_point()+ facet_wrap(~Target)+ theme_bw()+ 
  labs(x = "SREBF2 expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))


############################################################################################################################################################################################################
############################################################################################################################################################################################################
SCD <- geneexpression[geneexpression$Description=="SCD",]
SCD<-t(SCD)
str(SCD)
## SREBP1 is generated as a matrix somehow, and it needs to be converted into data frame
SCD <- as.data.frame(SCD)
setDT(SCD, keep.rownames = TRUE)[]
colnames(SCD) <- c("CCLE name","SCD")
#A one line option is: df$names<-rownames(df)
SCDexpression <- merge(SCD,Annotations,by='CCLE name')
## somehow the numbers of SCD columns are all changed into character 
SCDskin <- subset(SCDexpression, `Site Primary`=="skin")
write.csv(SCDskin,"SCDskin.csv")


SCDexpression_drug_data <- merge(SCDexpression,drug_data_short,by='CCLE name')
class(SCDexpression_drug_data$Compound)
as.factor(SCDexpression_drug_data$Compound)
class(SCDexpression_drug_data$SCD)
class(SCDexpression_drug_data$ActArea)
SCDexpression_drug_data$SCD <- as.numeric(as.character(SCDexpression_drug_data$SCD))
ggplot(SCDexpression_drug_data, aes(x=SCD, y=ActArea, colour=factor(Compound))) + geom_point()+ facet_wrap(~Compound)+ theme_bw()+ 
  labs(x = "SCD expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))

class(SCDexpression_drug_data$Compound)
SCDexpression_drug_data$Compound <- as.factor(SCDexpression_drug_data$Compound)
SCDexpression_drug_data_skin_PD <- subset(SCDexpression_drug_data, `Compound`=="Irinotecan")
cor.test(SCDexpression_drug_data_skin_PD$SCD, SCDexpression_drug_data_skin_PD$ActArea, method = "pearson")

SCDexpression_drug_data_17AAG <- subset(SCDexpression_drug_data, `Compound`=="17-AAG")
cor.test(SCDexpression_drug_data_17AAG$SCD, SCDexpression_drug_data_17AAG$ActArea, method = "pearson")

SCDexpression_drug_data_Paclitaxel <- subset(SCDexpression_drug_data, `Compound`=="Paclitaxel")
cor.test(SCDexpression_drug_data_Paclitaxel$SCD, SCDexpression_drug_data_Paclitaxel$ActArea, method = "pearson")

SCDexpression_drug_data_AZD6244 <- subset(SCDexpression_drug_data, `Compound`=="AZD6244")
cor.test(SCDexpression_drug_data_AZD6244$SCD, SCDexpression_drug_data_AZD6244$ActArea, method = "pearson")

SCDexpression_drug_data_PD0325901 <- subset(SCDexpression_drug_data, `Compound`=="PD-0325901")
cor.test(SCDexpression_drug_data_PD0325901$SCD, SCDexpression_drug_data_PD0325901$ActArea, method = "pearson")

SCDexpression_drug_data_MEK <- subset(SCDexpression_drug_data, `Target`=="MEK")
cor.test(SCDexpression_drug_data_MEK$SCD, SCDexpression_drug_data_MEK$ActArea, method = "pearson")

ggplot(SCDexpression_drug_data, aes(x=SCD, y=ActArea, colour=factor(Target))) + geom_point()+ facet_wrap(~Target)+ theme_bw()+ 
  labs(x = "SCD expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))


SCDexpression_drug_data_skin <- subset(SCDexpression_drug_data, `Site Primary`=="skin")
ggplot(SCDexpression_drug_data_skin, aes(x=SCD, y=ActArea, colour=factor(Compound))) + geom_point()+ facet_wrap(~Compound)+ theme_bw()+ 
  labs(x = "SCD expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))

class(SCDexpression_drug_data_skin$Compound)
SCDexpression_drug_data_skin$Compound <- as.factor(SCDexpression_drug_data_skin$Compound)
SCDexpression_drug_data_skin_MEK1 <- subset(SCDexpression_drug_data_skin, `Compound`=="AZD6244")
cor.test(SCDexpression_drug_data_skin_MEK1$SCD, SCDexpression_drug_data_skin_MEK1$ActArea, method = "pearson")

SCDexpression_drug_data_skin_MEK2 <- subset(SCDexpression_drug_data_skin, `Compound`=="PD-0325901")
cor.test(SCDexpression_drug_data_skin_MEK2$SCD, SCDexpression_drug_data_skin_MEK2$ActArea, method = "pearson")

SCDexpression_drug_data_skin_MEK <- subset(SCDexpression_drug_data_skin, `Target`=="MEK")
cor.test(SCDexpression_drug_data_skin_MEK$SCD, SCDexpression_drug_data_skin_MEK$ActArea, method = "pearson")

SCDexpression_drug_data_skin_HSP90 <- subset(SCDexpression_drug_data_skin, `Compound`=="17-AAG")
cor.test(SCDexpression_drug_data_skin_HSP90$SCD, SCDexpression_drug_data_skin_HSP90$ActArea, method = "pearson")

SCDexpression_drug_data_Topotecan <- subset(SCDexpression_drug_data, `Compound`=="Topotecan")
cor.test(SCDexpression_drug_data_Topotecan$SCD, SCDexpression_drug_data_Topotecan$ActArea, method = "pearson")

SCDexpression_drug_data_Paclitaxel <- subset(SCDexpression_drug_data, `Compound`=="Paclitaxel")
cor.test(SCDexpression_drug_data_Paclitaxel$SCD, SCDexpression_drug_data_Paclitaxel$ActArea, method = "pearson")

ggplot(SCDexpression_drug_data_skin, aes(x=SCD, y=ActArea, colour=factor(Target))) + geom_point()+ facet_wrap(~Target)+ theme_bw()+ 
  labs(x = "SCD expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))


############################################################################################################################################################################################################
############################################################################################################################################################################################################
FASN <- geneexpression[geneexpression$Description=="FASN",]
FASN<-t(FASN)
str(FASN)
## SREBP1 is generated as a matrix somehow, and it needs to be converted into data frame
FASN <- as.data.frame(FASN)
setDT(FASN, keep.rownames = TRUE)[]
colnames(FASN) <- c("CCLE name","FASN")
#A one line option is: df$names<-rownames(df)
FASNexpression <- merge(FASN,Annotations,by='CCLE name')
## somehow the numbers of SREBF1 columns are all changed into character 
FASNskin <- subset(FASNexpression, `Site Primary`=="skin")
write.csv(FASNskin,"FASNskin.csv")

FASNexpression_drug_data <- merge(FASNexpression,drug_data_short,by='CCLE name')
class(FASNexpression_drug_data$Compound)
as.factor(FASNexpression_drug_data$Compound)
class(FASNexpression_drug_data$SCD)
class(FASNexpression_drug_data$ActArea)
FASNexpression_drug_data$FASN <- as.numeric(as.character(FASNexpression_drug_data$FASN))
FASNexpression_drug_data$ActArea <- as.numeric(as.character(FASNexpression_drug_data$ActArea))

ggplot(FASNexpression_drug_data, aes(x=FASN, y=ActArea, colour=factor(Compound))) + geom_point()+ facet_wrap(~Compound)+ theme_bw()+ 
  labs(x = "FASN expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))

ggplot(FASNexpression_drug_data, aes(x=FASN, y=ActArea, colour=factor(Target))) + geom_point()+ facet_wrap(~Target)+ theme_bw()+ 
  labs(x = "FASN expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))

FASNexpression_drug_data_17AAG <- subset(FASNexpression_drug_data, `Compound`=="17-AAG")
cor.test(FASNexpression_drug_data_17AAG$FASN, FASNexpression_drug_data_17AAG$ActArea, method = "pearson")

FASNexpression_drug_data_Paclitaxel <- subset(FASNexpression_drug_data, `Compound`=="Paclitaxel")
cor.test(FASNexpression_drug_data_Paclitaxel$FASN, FASNexpression_drug_data_Paclitaxel$ActArea, method = "pearson")

FASNexpression_drug_data_Topotecan <- subset(FASNexpression_drug_data, `Compound`=="Topotecan")
cor.test(FASNexpression_drug_data_Topotecan$FASN, FASNexpression_drug_data_Topotecan$ActArea, method = "pearson")

FASNexpression_drug_data_AZD6244 <- subset(FASNexpression_drug_data, `Compound`=="AZD6244")
cor.test(FASNexpression_drug_data_AZD6244$FASN, FASNexpression_drug_data_AZD6244$ActArea, method = "pearson")

FASNexpression_drug_data_PD0325901 <- subset(FASNexpression_drug_data, `Compound`=="PD-0325901")
cor.test(FASNexpression_drug_data_PD0325901$FASN, FASNexpression_drug_data_PD0325901$ActArea, method = "pearson")

FASNexpression_drug_data_MEK <- subset(FASNexpression_drug_data, `Target`=="MEK")
cor.test(FASNexpression_drug_data_MEK$FASN, FASNexpression_drug_data_MEK$ActArea, method = "pearson")


FASNexpression_drug_data_skin <- subset(FASNexpression_drug_data, `Site Primary`=="skin")
ggplot(FASNexpression_drug_data_skin, aes(x=FASN, y=ActArea, colour=factor(Compound))) + geom_point()+ facet_wrap(~Compound)+ theme_bw()+ 
  labs(x = "FASN expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))

ggplot(FASNexpression_drug_data_skin, aes(x=FASN, y=ActArea, colour=factor(Target))) + geom_point()+ facet_wrap(~Target)+ theme_bw()+ 
  labs(x = "FASN expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))

class(FASNexpression_drug_data_skin$Compound)
FASNexpression_drug_data_skin$Compound <- as.factor(FASNexpression_drug_data_skin$Compound)
FASNexpression_drug_data_skin_AZD6244 <- subset(FASNexpression_drug_data_skin, `Compound`=="AZD6244")
cor.test(FASNexpression_drug_data_skin_AZD6244$FASN, FASNexpression_drug_data_skin_AZD6244$ActArea, method = "pearson")

FASNexpression_drug_data_skin_PD0325901 <- subset(FASNexpression_drug_data_skin, `Compound`=="PD-0325901")
cor.test(FASNexpression_drug_data_skin_PD0325901$FASN, FASNexpression_drug_data_skin_PD0325901$ActArea, method = "pearson")

FASNexpression_drug_data_skin_MEK <- subset(FASNexpression_drug_data_skin, `Target`=="MEK")
cor.test(FASNexpression_drug_data_skin_MEK$FASN, FASNexpression_drug_data_skin_MEK$ActArea, method = "pearson")

############################################################################################################################################################################################################
############################################################################################################################################################################################################
MITF <- geneexpression[geneexpression$Description=="MITF",]
MITF<-t(MITF)
str(MITF)
## SREBP1 is generated as a matrix somehow, and it needs to be converted into data frame
MITF <- as.data.frame(MITF)
setDT(MITF, keep.rownames = TRUE)[]
colnames(MITF) <- c("CCLE name","MITF")
#A one line option is: df$names<-rownames(df)
MITFexpression <- merge(MITF,Annotations,by='CCLE name')
## somehow the numbers of SCD columns are all changed into character 
MITFskin <- subset(MITFexpression, `Site Primary`=="skin")
write.csv(MITFskin,"MITFskin.csv")

MITFexpression_drug_data <- merge(MITFexpression,drug_data_short,by='CCLE name')
class(MITFexpression_drug_data$Compound)
as.factor(MITFexpression_drug_data$Compound)
class(MITFexpression_drug_data$MITF)
class(MITFexpression_drug_data$ActArea)
MITFexpression_drug_data$MITF <- as.numeric(as.character(MITFexpression_drug_data$MITF))
ggplot(MITFexpression_drug_data, aes(x=MITF, y=ActArea, colour=factor(Compound))) + geom_point()+ facet_wrap(~Compound)+ theme_bw()+ 
  labs(x = "MITF expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))

ggplot(MITFexpression_drug_data, aes(x=MITF, y=ActArea, colour=factor(Target))) + geom_point()+ facet_wrap(~Target)+ theme_bw()+ 
  labs(x = "MITF expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))


MITFexpression_drug_data_skin <- subset(MITFexpression_drug_data, `Site Primary`=="skin")
ggplot(MITFexpression_drug_data_skin, aes(x=MITF, y=ActArea, colour=factor(Compound))) + geom_point()+ facet_wrap(~Compound)+ theme_bw()+ 
  labs(x = "MITF expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))


ggplot(MITFexpression_drug_data_skin, aes(x=MITF, y=ActArea, colour=factor(Target))) + geom_point()+ facet_wrap(~Target)+ theme_bw()+ 
  labs(x = "MITF expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))


############################################################################################################################################################################################################
############################################################################################################################################################################################################
BRAF <- geneexpression[geneexpression$Description=="BRAF",]
EGFR<-t(EGFR)
str(EGFR)
## SREBP1 is generated as a matrix somehow, and it needs to be converted into data frame
EGFR <- as.data.frame(EGFR)
setDT(EGFR, keep.rownames = TRUE)[]
colnames(EGFR) <- c("CCLE name","EGFR")
#A one line option is: df$names<-rownames(df)
EGFRexpression <- merge(EGFR,Annotations,by='CCLE name')
## somehow the numbers of SCD columns are all changed into character 
EGFRskin <- subset(EGFRexpression, `Site Primary`=="skin")
write.csv(EGFRskin,"EGFRskin.csv")

EGFRexpression_drug_data <- merge(EGFRexpression,drug_data_short,by='CCLE name')
class(EGFRexpression_drug_data$Compound)
as.factor(EGFRexpression_drug_data$Compound)
class(EGFRexpression_drug_data$EGFR)
class(EGFRexpression_drug_data$ActArea)
EGFRexpression_drug_data$EGFR <- as.numeric(as.character(EGFRexpression_drug_data$EGFR))
ggplot(EGFRexpression_drug_data, aes(x=EGFR, y=ActArea, colour=factor(Compound))) + geom_point()+ facet_wrap(~Compound)+ theme_bw()+ 
  labs(x = "EGFR expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))


ggplot(SOX10expression_drug_data, aes(x=SOX10, y=ActArea, colour=factor(Target))) + geom_point()+ facet_wrap(~Target)+ theme_bw()+ 
  labs(x = "SOX10 expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))


SOX10expression_drug_data_skin <- subset(SOX10expression_drug_data, `Site Primary`=="skin")
ggplot(SOX10expression_drug_data_skin, aes(x=SOX10, y=ActArea, colour=factor(Compound))) + geom_point()+ facet_wrap(~Compound)+ theme_bw()+ 
  labs(x = "SOX10 expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))


ggplot(EGFRexpression_drug_data_skin, aes(x=EGFR, y=ActArea, colour=factor(Target))) + geom_point()+ facet_wrap(~Target)+ theme_bw()+ 
  labs(x = "EGFR expression log2(RMA)", y = "Activity Area") +   
  theme(axis.title=element_text(face="bold",size=12,color="black"),axis.text=element_text(size=12,face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))



############################################################################################################################################################################################################
############################################################################################################################################################################################################

## Merging more than 2 dataframes in R by rownames "CCLE name"
library(plyr)
SREBFSCDFASN <- join_all(list(SREBF1skin,SREBF2skin,SCDskin,FASNskin), by = 'CCLE name', type = 'full')
SREBFFASNSCDexpression <- subset(SREBFSCDFASN, select=c('Cell line primary name','SREBF1','SREBF2','SCD','FASN'))
write.csv(SREBFFASNSCDexpression,"SREBFFASNSCDexpression.csv")


## Tidying and normalising the data
SREBF1FASNSCD <- subset(SREBF1FASNSCDexpression, select=c('SCD','FASN','SREBF1'))
row.names(SREBF1FASNSCD) <- SREBF1FASNSCDexpression$'CCLE name'
class(SREBF1FASNSCD$SCD)
SREBF1FASNSCD$SCD <- as.numeric(as.character(SREBF1FASNSCD$SCD))
SREBF1FASNSCD$FASN <- as.numeric(as.character(SREBF1FASNSCD$FASN))
SREBF1FASNSCD$SREBF1<- as.numeric(as.character(SREBF1FASNSCD$SREBF1))
SREBF1FASNSCD_matrix <- as.matrix(SREBF1FASNSCD)
row.names(SREBF1FASNSCD_matrix) <- SREBF1FASNSCDexpression$'CCLE name'

# Normalize tidy data, transpose for row wise normalisation
dat.3 <- scale(SREBF1FASNSCD_matrix)

# Put data back in original form
dat.t3 <- t(dat.3)

# Check means are zero and std devs are 1
round(colMeans(dat.3),1)
apply(dat.3,2,sd)

## Hierarchical clustering
## Relating data points as distances
## Calculate distance between cell lines in rows
d1 <- dist(dat.3,method = "euclidean", diag = FALSE, upper = FALSE)
round(d1,3)
## Calculate distance between genes in column 
d2 <- dist(dat.t3,method = "euclidean", diag = FALSE, upper = TRUE)

round(d2,3)
# Clustering distance between experiments using Ward linkage
c1 <- hclust(d1, method = "ward.D2", members = NULL)
# Clustering distance between proteins using Ward linkage
c2 <- hclust(d2, method = "ward.D2", members = NULL)

# Check clustering by plotting dendrograms
par(mfrow=c(2,1),cex=0.5) # Make 2 rows, 1 col plot frame and shrink labels
plot(c1); plot(c2) # Plot both cluster dendrograms

# Set colours for heatmap, 25 increments
my_palette <- colorRampPalette(c("blue","white","red"))(n = 25)

# Plot heatmap with heatmap.2
par(cex.main=0.75) # Shrink title fonts on plot
heatmap.2(dat.3,                     # Tidy, normalised data
         # Colv=as.dendrogram(c2),     # Experiments clusters in cols
          Rowv=as.dendrogram(c1),     # Protein clusters in rows
          density.info="histogram",   # Plot histogram of data and colour key
          trace="none",               # Turn of trace lines from heat map
          col = my_palette,           # Use my colour scheme
          cexRow=0.5,cexCol=0.75)     # Amend row and column label fonts



heatmap(SREBF1FASNSCD_matrix, Colv=F, scale='none')
hc.rows <- hclust(dist(SREBF1FASNSCD_matrix))
plot(hc.rows)
hc.cols <- hclust(dist(t(SREBF1FASNSCD_matrix)))
heatmap(SREBF1FASNSCD_matrix[cutree(hc.rows,k=2)==1,], Colv=as.dendrogram(hc.cols), scale='none')
heatmap(SREBF1FASNSCD_matrix[cutree(hc.rows,k=2)==2,], Colv=as.dendrogram(hc.cols), scale='none')

library(RColorBrewer)
heatmap.2(SREBF1FASNSCD_matrix, trace="none", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

library(gtools)
distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method="average")
heatmap.2(SREBF1FASNSCD_matrix, trace="none", scale="row", zlim=c(-3,3), 
          col=cols, symbreak=FALSE) 
