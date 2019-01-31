library(data.table)
library(stringr)
library(ggplot2)
setwd("~/GDSC")


## http://www.cancerrxgene.org/downloads
## data from 'A landscape of pharmacogenomic interactions in cancer', Iorio F et al. Cell. 2016
GDSC.expression <- fread("Cell_line_RMA_proc_basalExp.txt",showProgress = T)
GDSC.annotation <- fread("Cell_Lines_Annotation.csv", showProgress = T)
GDSC.drugrespose <- fread("v17_fitted_dose_response.csv", showProgress = T)
GDSC.compounds <- fread("Screened_Compounds.csv", showProgress = T)

################################################################################
################################################################################
## boxplot to show individual gene expression across all cancer cell lines 
## (grouped by their tissue types)
plotGene <- function(gene.name) { 
  gene.object <- GDSC.expression[GDSC.expression$GENE_SYMBOLS == gene.name,]
  gene.object <- t(gene.object)
  class(gene.object)
  str(gene.object)
  ## gene.object is generated as a matrix somehow, and it needs to be converted into data frame
  gene.object <- as.data.frame(gene.object)
  setDT(gene.object, keep.rownames = TRUE)[]
  colnames(gene.object) <- c("COSMIC identifier","gene.name")
  gene.object <- gene.object[-c(1,2),]
  gene.object$'COSMIC identifier' <- str_replace_all(gene.object$'COSMIC identifier', "DATA.", "")
  GDSC.annotation$`COSMIC identifier` <- as.character(GDSC.annotation$`COSMIC identifier`)
  ## COSMIC identifier in GDSC.anotation file lists all cell line name as integer, 
  ## need to transform them into character before merge with expression data
  gene.object.expression <- merge(gene.object, GDSC.annotation, 
                                  by = 'COSMIC identifier')
  str(gene.object.expression)
# now SREBF1 expression values are classed as factors rather than numeric. 
# use the following script to convert factor into numeric safely
  gene.object.expression$gene.name <- as.numeric(levels(gene.object.expression$gene.name))[gene.object.expression$gene.name]
  gene.object.expression <- gene.object.expression[, -c(4,5,6,7,8)]
  gene.object.expression <- gene.object.expression[, c(1,2,3,4,5)]
  names(gene.object.expression) <- c("COSMIC_ID", "gene.name", 
                                     "Sample.Name", "Tissue.descriptor.1",
                                     "Tissue.descriptor.2")
  gene.object.expression <- as.data.frame(gene.object.expression)
  median1 <- within(gene.object.expression, 
                    Tissue.descriptor.1 <-  reorder(Tissue.descriptor.1, 
                                                    gene.name, median))
  ggplot(median1, aes(x = Tissue.descriptor.1, 
                      y = gene.name)) + 
    geom_boxplot()+ 
    theme_bw()+ 
    labs(x = "Primary Site (GDSC)",
         y = paste("RMA-normalized", gene.name, "mRNA expression")) +
    theme(axis.title = element_text(face = "bold",
                                    size = 12,
                                    color = "black"),
          axis.text = element_text(size  = 12,
                                   angle = 90, 
                                   hjust = 1, 
                                   face  = "bold",
                                   color = "black"),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"),
          panel.grid  = element_blank(),
          strip.text  = element_text(face  = "bold", 
                                     size  = 12, 
                                     color = "black"),
          legend.text = element_text(size  = 12, 
                                     face  = "bold",
                                     color = 'black'),
          legend.title = element_text(face = "bold", 
                                      size = 12, 
                                      color = "black"),
          legend.justification = c(1,1)) + 
    guides(fill = guide_legend(title = NULL))
  }

plotGene(gene.name = "EIF4E")
plotGene(gene.name = "EIF4G2")
plotGene(gene.name = "EIF4EBP1")
plotGene(gene.name = "RPS6KB1")
plotGene(gene.name = "MYC")


################################################################################
################################################################################
## generate dotplot overlaid on boxplot to show individual gene expression across all cancer cell line (grouped by their tissue types)
plotGene <- function(gene.name) { 
  gene.object <- GDSC.expression[GDSC.expression$GENE_SYMBOLS==gene.name,]
  gene.object<-t(gene.object)
  class(gene.object)
  str(gene.object)
  ## gene.object is generated as a matrix somehow, and it needs to be converted into data frame
  gene.object <- as.data.frame(gene.object)
  setDT(gene.object, keep.rownames = TRUE)[]
  colnames(gene.object) <- c("COSMIC identifier","gene.name")
  gene.object <- gene.object[-c(1,2),]
  gene.object$'COSMIC identifier' <- str_replace_all(gene.object$'COSMIC identifier', "DATA.", "")
  GDSC.annotation$`COSMIC identifier` <- as.character(GDSC.annotation$`COSMIC identifier`)
  ## COSMIC identifier in GDSC.anotation file lists all cell line name as integer, 
  ## need to transform them into character before merge with expression data
  gene.object.expression <- merge(gene.object, GDSC.annotation,by='COSMIC identifier')
  str(gene.object.expression)
  # now SREBF1 expression values are classed as factors rather than numeric. 
  # use the following script to convert factor into numeric safely
  gene.object.expression$gene.name <- as.numeric(levels(gene.object.expression$gene.name))[gene.object.expression$gene.name]
  gene.object.expression <- gene.object.expression[, -c(4,5,6,7,8)]
  gene.object.expression <- gene.object.expression[, c(1,2,3,4,5)]
  names(gene.object.expression) <- c("COSMIC_ID", "gene.name", "Sample.Name", "Tissue.descriptor.1","Tissue.descriptor.2")
  gene.object.expression <- as.data.frame(gene.object.expression)
  median1 <- within(gene.object.expression, Tissue.descriptor.1 <-  reorder(Tissue.descriptor.1, gene.name, mean))
  ggplot(median1, aes(x=Tissue.descriptor.1, y=gene.name)) + 
    geom_boxplot(outlier.color = NA, width=0.4)+
    geom_dotplot(binaxis="y", binwidth=.05,stackdir="centerwhole")+ theme_bw()+ 
    labs(x = "Primary Site (GDSC)",
         y = paste("RMA-normalized", gene.name, "mRNA expression")) +
    theme(axis.title=element_text(face="bold",size=12,color="black"),
          axis.text=element_text(size=12,angle = 90, hjust = 1, face="bold",color="black"),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 12, colour = "black"),
          legend.text = element_text(size=12,face="bold",colour = 'black'),
          legend.title = element_text(face = "bold", size = 12, colour = "black"),
          legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))
}

plotGene(gene.name="FASN")
plotGene(gene.name="SCD")
plotGene(gene.name="ACLY")
plotGene(gene.name="SREBF1")
plotGene(gene.name="HMGCR")
plotGene(gene.name="EIF4EBP1")
plotGene(gene.name="RPS6KB1")
plotGene(gene.name="EIF4G1")
plotGene(gene.name="EIF4E")

## or use apply function 
x=c("FASN", "SCD","ACLY","EIF4A1","EIF4E","EIF4G1","EIF4EBP1","RPS6KB1","MYC")
lapply(x, plotGene)

############################################################################################################################################################################################################
############################################################################################################################################################################################################
x=c("FASN", "SCD","ACLY","SREBF1","HMGCR","HMGCS1","LDLR","SREBF2","MITF")
gene.object <- NULL
for (i in x) gene.object <- rbind(gene.object, GDSC.expression[GDSC.expression$GENE_SYMBOL== i,])

str(gene.object)
gene.object <- t(gene.object)
gene.object <- as.data.frame(gene.object)
# change gene.object into data.frame is important for rowname change!
setDT(gene.object, keep.rownames = TRUE)[]
colnames(gene.object) <- c("COSMIC identifier",x)
gene.object <- gene.object[-c(1,2),]
gene.object$'COSMIC identifier' <- str_replace_all(gene.object$'COSMIC identifier', "DATA.", "")
GDSC.annotation$`COSMIC identifier` <- as.character(GDSC.annotation$`COSMIC identifier`)
## COSMIC identifier in GDSC.anotation file lists all cell line name as integer, 
## need to transform them into character before merge with expression data
gene.object.expression <- merge(gene.object, GDSC.annotation,by='COSMIC identifier')
## $i will not evaluate i to find out it's current value. 
## Should instead access and assign to the column using double brackets, 
## double brackets [[]] work when i is a name or index number:
for (i in x) {
  gene.object.expression[[i]] <- as.numeric(levels(gene.object.expression[[i]]))[gene.object.expression[[i]]]
}
gene.object.expression <- gene.object.expression[, -c(12:16)]
gene.object.expression <- gene.object.expression[, -c(14:17)]
names(gene.object.expression) <- c("COSMIC_ID", x, "Sample.Name", "Tissue.descriptor.1","Tissue.descriptor.2")
gene.object.expression <- as.data.frame(gene.object.expression)
gene.object.expression$Tissue.descriptor.1 <- as.factor(gene.object.expression$Tissue.descriptor.1)

median <- within(gene.object.expression, GDSCTissue.descriptor.1 <-  reorder(GDSCTissue.descriptor.1, SCD, median))

p <- ggplot(gene.object.expression, aes(FASN, SCD))  
p + geom_point()
p + geom_point(aes(colour= factor(Tissue.descriptor.1)))+ facet_wrap(~Tissue.descriptor.1)  

p <- ggplot(gene.object.expression, aes(MITF, SCD))  
p + geom_point()
p + geom_point(aes(colour= factor(Tissue.descriptor.1)))+ facet_wrap(~Tissue.descriptor.1) 


cor.test(gene.object.expression$FASN,gene.object.expression$SCD, method = "pearson" )                   
cor.test(gene.object.expression$SREBF1,gene.object.expression$SCD, method = "pearson" )                   
cor.test(gene.object.expression$SREBF2,gene.object.expression$SCD, method = "pearson" ) 
cor.test(gene.object.expression$HMGCR,gene.object.expression$SCD, method = "pearson" ) 
cor.test(gene.object.expression$SCD,gene.object.expression$MITF, method = "pearson" ) 

############################################################################################################################################################################################################
## combine  gene expression with drug sensitivity data
head(GDSC.compounds)
names(GDSC.compounds) <- c("DRUG_ID", "DRUG_NAME", "SYNONYMS", "TARGET","TARGET_PATHWAY")
GDSC.drugrespose.annotation <- merge(GDSC.drugrespose, GDSC.compounds,by='DRUG_ID')
GDSC.drugrespose.annotation$TARGET_PATHWAY <- as.factor(GDSC.drugrespose.annotation$TARGET_PATHWAY)  

plotDrug <- function(gene.name) { 
  ## gene.name "gene of interest"
  gene.object <- GDSC.expression[GDSC.expression$GENE_SYMBOLS==gene.name,]
  gene.object<-t(gene.object)
  class(gene.object)
  str(gene.object)
  ## gene.object is generated as a matrix somehow, and it needs to be converted into data frame
  gene.object <- as.data.frame(gene.object)
  setDT(gene.object, keep.rownames = TRUE)[]
  colnames(gene.object) <- c("COSMIC identifier", gene.name)
  gene.object <- gene.object[-c(1,2),]
  gene.object$'COSMIC identifier' <- str_replace_all(gene.object$'COSMIC identifier', "DATA.", "")
  GDSC.annotation$`COSMIC identifier` <- as.character(GDSC.annotation$`COSMIC identifier`)
  ## COSMIC identifier in GDSC.annotation file lists all cell line name as integer, 
  ## need to transform them into character before merge with expression data
  gene.object <- merge(gene.object, GDSC.annotation,by='COSMIC identifier')
  gene.object <- gene.object[, -c(4,5,6,7,8)]
  gene.object <- gene.object[, c(1,2,3,4,5)]
  str(gene.object)
  # now SREBF1 expression values are classed as factors rather than numeric. 
  # use the following script to convert factor into numeric safely
  gene.object$gene.name <- as.numeric(levels(gene.object$gene.name))[gene.object$gene.name]
  names(gene.object) <- c("COSMIC_ID", gene.name, "Sample.Name", "Tissue.descriptor.1","Tissue.descriptor.2")
  gene.object <- as.data.frame(gene.object) 
  
  expression_drugdata <- merge(gene.object.expression,GDSC.drugrespose.annotation,by='COSMIC_ID')
  expression_drugdata_skin <- subset(expression_drugdata, `Tissue.descriptor.1`=="skin")
  expression_drugdata_skin_goi <- subset(expression_drugdata_skin, `TARGET_PATHWAY`==gene.name)
  ggplot(expression_drugdata_skin_goi, aes(x=log2(SCD), y=log2(RMSE), colour=factor(TARGET_PATHWAY))) + geom_point()+ facet_wrap(~TARGET)+ theme_bw()+ 
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
x <- levels(GDSC.drugrespose.annotation$TARGET_PATHWAY)
lapply(x, plotDrug)
  
gene.name= "FASN"  
  
  
  expression_drugdata_ERK <- subset(expression_drugdata, `TARGET_PATHWAY`=="ERK MAPK signaling")
  expression_drugdata_skin <- subset(expression_drugdata, `Tissue.descriptor.1`=="skin")
  expression_drugdata_skin_ERK <- subset(expression_drugdata_skin, `TARGET_PATHWAY`=="ERK MAPK signaling")
  expression_drugdata_skin_JNK <- subset(expression_drugdata_skin, `TARGET_PATHWAY`=="JNK and p38 signaling")
  str(expression_drugdata$TARGET_PATHWAY)
  
  str(expression_drugdata)
  str(expression_drugdata_skin)
  expression_drugdata$TARGET <- as.factor(expression_drugdata$TARGET)
  expression_drugdata$TARGET_PATHWAY <- as.factor(expression_drugdata$TARGET_PATHWAY)
  expression_drugdata$DRUG_NAME <- as.factor(expression_drugdata$DRUG_NAME)
  goi <- levels(expression_drugdata$TARGET_PATHWAY)
}
############################################################################################################################################################################################################
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
GDCS.Cor <- t(cor(expression_drugdata_skin$"SCD", expression_drugdata_skin[,2:10]))
############################################################################################################################################################################################################
CorDrug <- function(goi) {
  expression_drugdata_skin_goi <- subset(expression_drugdata_skin, `DRUG_NAME`==goi)
  goi <- cor.test(expression_drugdata_skin_goi$SCD, expression_drugdata_skin_goi$RMSE, method = "pearson")
  data.frame(goi[c("estimate","p.value","statistic")])
}

class(expression_drugdata_skin$DRUG_NAME)
expression_drugdata_skin$DRUG_NAME <- as.factor(expression_drugdata_skin$DRUG_NAME)
x <- levels(expression_drugdata_skin$DRUG_NAME)
CorDrugResult <- t(sapply(x, CorDrug))
CorDrugResult <- as.data.frame(CorDrugResult)
str(CorDrugResult$p.value)
CorDrugResult$statistic <- as.numeric(CorDrugResult$statistic)
write.csv(CorDrugResult, "CorDrugResult.csv", fileEncoding = "utf8")

CorTarget <- function(goi) {
  expression_drugdata_skin_goi <- subset(expression_drugdata_skin, `TARGET`==goi)
  goi <- cor.test(expression_drugdata_skin_goi$SCD, expression_drugdata_skin_goi$RMSE, method = "pearson")
  data.frame(goi[c("estimate","p.value","statistic")])
}
class(expression_drugdata_skin$TARGET)
expression_drugdata_skin$TARGET <- as.factor(expression_drugdata_skin$TARGET)
y <- levels(expression_drugdata_skin$TARGET)
CorTARGETResult <- t(sapply(y, CorTarget))
CorTARGETResult <- as.data.frame(CorTARGETResult)
CorTARGETResult$estimate <- as.numeric(CorTARGETResult$estimate)
CorTARGETResult$statistic <- as.numeric(CorTARGETResult$statistic)

CorTarget_PATHWAY <- function(goi) {
  expression_drugdata_skin_goi <- subset(expression_drugdata_skin, `TARGET_PATHWAY`==goi)
  goi <- cor.test(expression_drugdata_skin_goi$SCD, expression_drugdata_skin_goi$RMSE, method = "pearson")
  data.frame(goi[c("estimate","p.value","statistic")])
}
z <- levels(expression_drugdata$TARGET_PATHWAY)
CorTARGET_PATHWAYResult <- t(sapply(z, CorTarget_PATHWAY))
CorTARGET_PATHWAYResult <- as.data.frame(CorTARGET_PATHWAYResult)
CorTARGET_PATHWAYResult$estimate <- as.numeric(CorTARGET_PATHWAYResult$estimate)
CorTARGET_PATHWAYResult$statistic <- as.numeric(CorTARGET_PATHWAYResult$statistic)