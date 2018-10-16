library(data.table)
library(ggplot2)

gene <- fread("GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.csv",showProgress = T)
Annotations <- fread("GTEx_Data_V6_Annotations_SampleAttributesDS.txt",showProgress = T)
Annotations <-Annotations[,c(1,7)]
## showProgress = T is necessary, otherwise "Error: isLOGICAL(showProgress) is not TRUE"

######################################################################################################
plotGene <- function(goi) {
  go <- gene[gene$Description==goi,]
  go<-t(go)
  ## SREBP1 is generated as a matrix somehow, and it needs to be converted into data frame
  go <- data.frame(go)
  setDT(go, keep.rownames = TRUE)[]
  colnames(go) <- c("SAMPID", "goi")
  #A one line option is: df$names<-rownames(df)
  goexpression <- merge(go,Annotations,by='SAMPID')
  ## somehow the numbers of SREBF1 columns are all changed into character 
  goexpression$SMTSD <- as.factor(goexpression$SMTSD)
  #  In particular, as.numeric applied to a factor is meaningless, and may happen by implicit coercion. 
  #  To transform a factor f to approximately its original numeric values, as.numeric(levels(f))[f] is recommended. 
  goexpression$goi <- as.numeric(levels(goexpression$goi))[goexpression$goi]
  goexpression <- as.data.frame(goexpression)
  ## draw boxplot for FASN expression across different tissues
  mean <- within(goexpression, SMTSD <-reorder(SMTSD, log2(goi), median))
  ggplot(mean, aes(x=SMTSD, y=log2(goi))) + geom_boxplot()+ theme_bw()+ 
    labs(x = "Tissue types (GTEx)",
         y = paste("log2(", goi, "RNA counts)")) +
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

plotGene(goi="FASN")
plotGene(goi="SCD")
plotGene(goi="SREBF1")
plotGene(goi="HMGCR")
plotGene(goi="HMGCS1")
plotGene(goi="SREBF2")
plotGene(goi="MITF")




