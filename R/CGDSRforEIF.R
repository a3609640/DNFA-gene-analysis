# the following script generates the plot and statistics
# for EIF expression and mutations data from TCGA provisional dataset.
# install.packages("cgdsr")
library(car)
library(cgdsr)
library(cowplot)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(httr)
library(plyr)
library(reshape2)
library(stringr)
library(survival)
library(survMisc)

# Create CGDS object
mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
# mycgds = CGDS("http://www.cbioportal.org/")
test(mycgds)
# Get list of cancer studies at server
getCancerStudies(mycgds)
# Get cases from TCGA provisional studies only
EIF.gene <- c("SCD","EIF4A1","EIF4E","EIF4G1","EIF4EBP1","RPS6KB1","MYC")
names(EIF.gene) <- EIF.gene

####################################################################
## plot EIF RNASeq data from TCGA provisional cancer study groups ##
####################################################################
plot.EIF.tcga <- function(EIF){
  # Get EIF RNAseq data from all TCGA study groups
  tcga.pro.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)", getCancerStudies(mycgds)$name), ]
  # "tcag_study_list" contains all the tcga cancer studies
  tcga.study.list <- tcga.pro.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  # use lappy to pull out all the caselists within tcga.study.list
  # because we named each elements in tcga.study.list,
  # lappy will return a large list, each element (with a cancer study name)
  # in that list is a data-table
  tcga.pro.caselist <- lapply(tcga.study.list, caselist)
  tcga.pro.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  # for example, tcga.pro.caselist[[1]] shows the dataframe of caselist
  # in laml study group.
  # to choose case_list_id that is labeled with laml_tcga_rna_seq_v2_mrna,
  # we use the following tcag_provisional_caselist[[1][8,1]
  # a <- tcga.pro.caselist[[1]][
  # grep("tcga_rna_seq_v2_mrna", tcga.pro.caselist[[1]]$case_list_id),
  # ][1,1]
  # b <- tcga.pro.geneticprofile[[1]][
  # grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
  # tcga.pro.geneticprofile[[1]]$genetic_profile_name), ][1,1]
  # how do we do this for all study groups from [[1]] to  [[32]]?
  caselist.RNAseq <- function(x) {
    tcga.pro.caselist[[x]][
      grep("tcga_rna_seq_v2_mrna",
           tcga.pro.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.RNAseq <- function(x) {
    tcga.pro.geneticprofile[[x]][
      # double backslash \\ suppress the special meaning of ( )
      # in regular expression
      grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
           tcga.pro.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
  # test the functions: caselist.RNAseq () and geneticprofile.RNAseq ()
  # caselist.RNAseq = caselist.RNAseq ('acc_tcga')
  # geneticprofile.RNAseq = geneticprofile.RNAseq ('acc_tcga')
  # Wrap two functions: geneticprofile.RNAseq(x), caselist.RNAseq(x)
  # within TCGA_ProfileData_RNAseq(x)
  tcga.profiledata.RNAseq <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
                   genename,
                   geneticprofile,
                   caselist)
  }
  EIF.tcga.RNAseq <- function(x, y) {
    tcga.profiledata.RNAseq(x,
                            geneticprofile.RNAseq(y),
                            caselist.RNAseq(y))
  }
  EIF.RNAseq.tcga.all <- function(x) {
    test <- lapply(tcga.study.list,
                   function(y) mapply(EIF.tcga.RNAseq, x, y))
    df2 <- melt(test)
    colnames(df2) <- c("RNAseq", "EIFgene", "TCGAstudy")
    df2 <- data.frame(df2)
  }
  df2 <- EIF.RNAseq.tcga.all(EIF)
  df2$EIFgene <- as.factor(df2$EIFgene)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  df2 <- na.omit(df2)
  # plot EIF gene expression across all TCGA groups ##
  m <- paste0(EIF, ".", EIF)
  mean <- within(df2[df2$EIFgene == m,], # TCGAstudy is one column in df2
                 TCGAstudy <- reorder(TCGAstudy, log2(RNAseq), median))
  print(
    ggplot(mean,
           aes(x     = TCGAstudy,
               y     = log2(RNAseq),
               color = TCGAstudy)) +
      geom_boxplot(alpha    = .01,
                   width    = .5,
                   position = position_dodge(width = .9)) +
      labs(x = "Tumor types (TCGA)",
           y = paste0("log2(", EIF, " RNA counts)")) +
      theme(axis.title  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            axis.text.x = element_text(size   = 9,
                                       angle  = 45,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.text.y = element_text(size   = 9,
                                       angle  = 0,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.line.x = element_line(color  = "black"),
            axis.line.y = element_line(color  = "black"),
            panel.grid  = element_blank(),
            strip.text  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            legend.position = "none"))
}

plot.EIF.tcga("EIF4E")
sapply(EIF.gene, plot.EIF.tcga)

#############################################################
## plot EIF RNASeq data from pan TCGA cancer study groups ##
#############################################################
plot.EIF.pan.tcga <- function(EIF){
  # Get EIF RNAseq data from all TCGA study groups
  tcga.pan.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
  # "tcag_study_list" contains all the tcga cancer studies
  tcga.study.list <- tcga.pan.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  # use lappy to pull out all the caselists within tcga.study.list
  # because we named each elements in tcga.study.list,
  # lappy will return a large list, each element (with a cancer study name)
  # in that list is a data-table
  tcga.pan.caselist <- lapply(tcga.study.list, caselist)
  tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  # for example, tcga.pro.caselist[[1]] shows the dataframe of caselist
  # in laml study group.
  # to choose case_list_id that is labeled with laml_tcga_rna_seq_v2_mrna,
  # we use the following tcag_provisional_caselist[[1][8,1]
  # a <- tcga.pro.caselist[[1]][
  # grep("tcga_rna_seq_v2_mrna", tcga.pro.caselist[[1]]$case_list_id),
  # ][1,1]
  # b <- tcga.pro.geneticprofile[[1]][
  # grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
  # tcga.pro.geneticprofile[[1]]$genetic_profile_name), ][1,1]
  # how do we do this for all study groups from [[1]] to  [[32]]?
  caselist.RNAseq <- function(x) {
    tcga.pan.caselist[[x]][
      grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
           tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.RNAseq <- function(x) {
    tcga.pan.geneticprofile[[x]][
      # double backslash \\ suppress the special meaning of ( )
      # in regular expression
      grep("mRNA Expression, RSEM",
           tcga.pan.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
  # test the functions: caselist.RNAseq () and geneticprofile.RNAseq ()
  # caselist.RNAseq = caselist.RNAseq ('acc_tcga')
  # geneticprofile.RNAseq = geneticprofile.RNAseq ('acc_tcga')
  # Wrap two functions: geneticprofile.RNAseq(x), caselist.RNAseq(x)
  # within TCGA_ProfileData_RNAseq(x)
  tcga.profiledata.RNAseq <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
                   genename,
                   geneticprofile,
                   caselist)
  }
  EIF.tcga.RNAseq <- function(x, y) {
    tcga.profiledata.RNAseq(x,
                            geneticprofile.RNAseq(y),
                            caselist.RNAseq(y))
  }
  EIF.RNAseq.tcga.all <- function(x) {
    test <- lapply(tcga.study.list,
                   function(y) mapply(EIF.tcga.RNAseq, x, y))
    df2 <- melt(test)
    colnames(df2) <- c("RNAseq", "EIFgene", "TCGAstudy")
    df2 <- data.frame(df2)
  }
  df2 <- EIF.RNAseq.tcga.all(EIF)
  df2$EIFgene <- as.factor(df2$EIFgene)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  df2 <- na.omit(df2)
  # plot EIF gene expression across all TCGA groups ##
  m <- paste0(EIF, ".", EIF)
  mean <- within(df2[df2$EIFgene == m,], # TCGAstudy is one column in df2
                 TCGAstudy <- reorder(TCGAstudy, log2(RNAseq), median))
  print(
    ggplot(mean,
           aes(x     = TCGAstudy,
               y     = log2(RNAseq),
               color = TCGAstudy)) +
      geom_boxplot(alpha    = .01,
                   width    = .5,
                   position = position_dodge(width = .9)) +
      labs(x = "Tumor types (TCGA)",
           y = paste0("log2(", EIF, " RNA counts)")) +
      theme(axis.title  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            axis.text.x = element_text(size   = 9,
                                       angle  = 45,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.text.y = element_text(size   = 9,
                                       angle  = 0,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.line.x = element_line(color  = "black"),
            axis.line.y = element_line(color  = "black"),
            panel.grid  = element_blank(),
            strip.text  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            legend.position = "none"))
}

plot.EIF.pan.tcga("EIF4E")
sapply(EIF.gene, plot.EIF.pan.tcga)

##########################################################
## plot RNAseq data of EIF complex in TCGA study groups ##
##########################################################
EIF.RNAseq.data <- getProfileData(mycgds,
                                  c("EIF4A1","EIF4B","EIF4E",
                                    "EIF4G1","EIF4EBP1","RPS6KB1","MYC"),
                                  "laml_tcga_rna_seq_v2_mrna",
                                  "laml_tcga_all")
EIF.RNAseq.data <- na.omit(EIF.RNAseq.data)
boxplot(log2(EIF.RNAseq.data), main="EIF RNAseq in laml")


EIF.RNAseq.data <- getProfileData(mycgds,
                                  c("EIF4A1","EIF4B","EIF4E",
                                    "EIF4G1","EIF4EBP1","RPS6KB1","MYC"),
                                  "acc_tcga_rna_seq_v2_mrna",
                                  "acc_tcga_all")
EIF.RNAseq.data <- na.omit(EIF.RNAseq.data)
boxplot(log2(EIF.RNAseq.data), main="EIF RNAseq data in hnsc")


EIF.RNAseq.data <- getProfileData(mycgds,
                                  c("EIF4A1","EIF4B","EIF4E",
                                    "EIF4G1","EIF4EBP1","RPS6KB1","MYC"),
                                  "skcm_tcga_rna_seq_v2_mrna",
                                  "skcm_tcga_all")
EIF.RNAseq.data <- na.omit(EIF.RNAseq.data)
boxplot(log2(EIF.RNAseq.data), main="EIF RNAseq data in skcm")


EIF.RNAseq.data <- getProfileData(mycgds,
                                  c("EIF4A1","EIF4B","EIF4E",
                                    "EIF4G1","EIF4EBP1","RPS6KB1","MYC"),
                                  "dlbc_tcga_rna_seq_v2_mrna",
                                  "dlbc_tcga_all")
EIF.RNAseq.data <- na.omit(EIF.RNAseq.data)
boxplot(log2(EIF.RNAseq.data), main="EIF RNAseq data in dlbc")

EIF.RNAseq.data <- getProfileData(mycgds,
                                  c("EIF4A1","EIF4B","EIF4E",
                                    "EIF4G1","EIF4EBP1","RPS6KB1","MYC"),
                                  "esca_tcga_rna_seq_v2_mrna",
                                  "esca_tcga_all")
EIF.RNAseq.data <- na.omit(EIF.RNAseq.data)
boxplot(log2(EIF.RNAseq.data), main="EIF RNAseq data in Esophageal cancer")

##########################################################
tcga.pan.studies <- getCancerStudies(mycgds)[
  grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
# "tcag_study_list" contains all the tcga cancer studies
tcga.study.list <- tcga.pan.studies$cancer_study_id
names(tcga.study.list) <- tcga.study.list
caselist <- function(x) getCaseLists(mycgds, x)
geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
# use lappy to pull out all the caselists within tcga.study.list
# because we named each elements in tcga.study.list,
# lappy will return a large list, each element (with a cancer study name)
# in that list is a data-table
tcga.pan.caselist <- lapply(tcga.study.list, caselist)
tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
# for example, tcga.pro.caselist[[1]] shows the dataframe of caselist
# in laml study group.
# to choose case_list_id that is labeled with laml_tcga_rna_seq_v2_mrna,
# we use the following tcag_provisional_caselist[[1][8,1]
# a <- tcga.pro.caselist[[1]][
# grep("tcga_rna_seq_v2_mrna", tcga.pro.caselist[[1]]$case_list_id),
# ][1,1]
# b <- tcga.pro.geneticprofile[[1]][
# grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
# tcga.pro.geneticprofile[[1]]$genetic_profile_name), ][1,1]
# how do we do this for all study groups from [[1]] to  [[32]]?
caselist.RNAseq <- function(x) {
  tcga.pan.caselist[[x]][
    grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
         tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
}
geneticprofile.RNAseq <- function(x) {
  tcga.pan.geneticprofile[[x]][
    # double backslash \\ suppress the special meaning of ( )
    # in regular expression
    grep("mRNA Expression, RSEM",
         tcga.pan.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
}

tcga.profiledata.RNAseq <- function(genename,
                                    geneticprofile,
                                    caselist) {
  getProfileData(mycgds,
                 genename,
                 geneticprofile,
                 caselist)
}

get.EIF.RNAseq.tcga <- function(x) {
  EIF.gene <- c("EIF4E","EIF4G1","EIF4EBP1","RPS6KB1")
  EIF.RNAseq.tcga <- tcga.profiledata.RNAseq(EIF.gene,
                                             geneticprofile.RNAseq(x),
                                             caselist.RNAseq(x))
  EIF.RNAseq.tcga$rn <- rownames(EIF.RNAseq.tcga)
  return(EIF.RNAseq.tcga)
}

get.EIF.score.tcga <- function(x){
  EIF.RNAseq.tcga <- get.EIF.RNAseq.tcga(x) 
  EIF.score.tcga <- EIF.RNAseq.tcga
  EIF.score.tcga$EIF4Escore <- EIF.RNAseq.tcga$EIF4E/EIF.RNAseq.tcga$EIF4E
  EIF.score.tcga$EIF4G1score <- EIF.RNAseq.tcga$EIF4G1/EIF.RNAseq.tcga$EIF4E
  EIF.score.tcga$EIF4EBP1score <- EIF.RNAseq.tcga$EIF4EBP1/EIF.RNAseq.tcga$EIF4E
#  EIF.score.tcga$RPS6KB1score <- EIF.RNAseq.tcga$RPS6KB1/EIF.RNAseq.tcga$EIF4E
  EIF.score.tcga <- EIF.score.tcga [, 5:8]
  return(EIF.score.tcga)
}
# EIF.score.tcga$GEBPscore <- EIF.RNAseq.tcga$EIF4G1/EIF.RNAseq.tcga$EIF4EBP1

plot.EIF.RNAseq.score <- function (x) {
  EIF.RNAseq.tcga <- get.EIF.RNAseq.tcga(x)
  EIF.score.tcga <- get.EIF.score.tcga(x)
  par(mfrow=c(1,2))
  boxplot(log2(EIF.RNAseq.tcga[, 
                               c("EIF4E", "EIF4G1", 
                                 "EIF4EBP1", "RPS6KB1")]),
          main= paste0("EIF RNAseq counts in ", x),
          las = 2)
  boxplot(log2(EIF.score.tcga[,
                              c("EIF4Escore","EIF4G1score",
                                "EIF4EBP1score","RPS6KB1score")]),
          main= paste0("EIF scores in ", x),
          las = 2)
}

lapply(tcga.study.list, plot.EIF.RNAseq.score)


### to be tested!
plot.EIF.score.all.tcga <- function(x) {
  EIF.score.tcga <- lapply(tcga.study.list, get.EIF.score.tcga)
  EIF.score.tcga.all.tumors <- melt(EIF.score.tcga)
  colnames(EIF.score.tcga.all.tumors) <- c("rn", "EIFgene", "Score", "TCGAstudy")
  EIF.score.tcga.all.tumors <- data.frame(EIF.score.tcga.all.tumors)
  EIF.score.tcga.all.tumors$EIFgene <- as.factor(EIF.score.tcga.all.tumors$EIFgene)
  EIF.score.tcga.all.tumors$TCGAstudy <- as.factor(EIF.score.tcga.all.tumors$TCGAstudy)
  ### to be tested!  
  #  EIF.score.tcga.all.tumors <- EIF.score.tcga.all.tumors[EIF.score.tcga.all.tumors$EIFgene == x,]
  #  mean <- within(EIF.score.tcga.all.tumors, TCGAstudy <- reorder(TCGAstudy, log2(x), median))
  y <- paste0(x, "score")
  median <- within(EIF.score.tcga.all.tumors[EIF.score.tcga.all.tumors$EIFgene == y,], # TCGAstudy is one column in df2
                   TCGAstudy <- reorder(TCGAstudy, log2(Score), median))
  print(
    ggplot(median,
           aes(x     = TCGAstudy,
               y     = log2(Score),
               color = TCGAstudy)) +
      geom_boxplot(alpha    = .01,
                   width    = .5,
                   position = position_dodge(width = .9)) +
      labs(x = "Tumor types (TCGA)",
           y = paste0("log2(EIF4G1/EIF4EBP1 ratio)")) +
      theme(axis.title  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            axis.text.x = element_text(size   = 9,
                                       angle  = 45,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.text.y = element_text(size   = 9,
                                       angle  = 0,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.line.x = element_line(color  = "black"),
            axis.line.y = element_line(color  = "black"),
            panel.grid  = element_blank(),
            strip.text  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            legend.position = "none"))
  
}


plot.EIF.score.all.tcga("GEBP")


##########################################################
##########################################################
plotEIF <-  function (x) {
  #  name <- deparse(substitute(x))
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
  ggplot(x,
         aes(x = variable,
             y = log2(value))) +
    geom_boxplot(alpha    = .01, 
                 size     = .75,
                 width    = .75,
                 position = position_dodge(width = .9)) +
    #    labs(title = paste0(name," n = 8555"),
    #         x     = "eIF4F complex components",
    #         y     = paste0("log2(value)")) +
    theme_bw() +
    theme(plot.title  = black_bold_tahoma_12,
          axis.title  = black_bold_tahoma_12,
          axis.text.x = black_bold_tahoma_12_45,
          axis.text.y = black_bold_tahoma_12,
          axis.line.x = element_line(color  = "black"),
          axis.line.y = element_line(color  = "black"),
          panel.grid  = element_blank(),
          strip.text  = black_bold_tahoma_12,
          legend.position = "none")
}

plot.EIFandScore.all.tumors <- function (){
  EIF.RNAseq.tcga <- lapply(tcga.study.list, get.EIF.RNAseq.tcga)
  EIF.score.tcga <- lapply(tcga.study.list, get.EIF.score.tcga)
  EIF.RNAseq.tcga.all.tumors <- melt(EIF.RNAseq.tcga)
  #  x1  = factor(x, levels=c("B", "C", "A"))
  #  levels(EIF.RNAseq.tcga.all.tumors$variable)
  EIF.RNAseq.tcga.all.tumors$variable <- ordered(EIF.RNAseq.tcga.all.tumors$variable, 
                                                 levels = c("EIF4E",
                                                            "EIF4G1",
                                                            "EIF4EBP1",
                                                            "RPS6KB1"))
  number <- nrow(EIF.RNAseq.tcga.all.tumors)/4
  EIF.score.tcga.all.tumors <- melt(EIF.score.tcga)
  my_comparison1 <- list( c("EIF4E", "EIF4G1"), 
                          c("EIF4G1", "EIF4EBP1"), 
                          c("EIF4E", "EIF4EBP1"),
                          c("EIF4E", "RPS6KB1"),
                          c("EIF4EBP1", "RPS6KB1"))
  my_comparison2 <- list( c("EIF4Escore", "EIF4G1score"), 
                          c("EIF4G1score", "EIF4EBP1score"), 
                          c("EIF4Escore", "EIF4EBP1score"),
                          c("EIF4Escore", "RPS6KB1score"),
                          c("EIF4EBP1score", "RPS6KB1score"))
  p1 <- plotEIF(EIF.RNAseq.tcga.all.tumors) +
    labs(title = paste0("All tumors n = ", number),
         x     = "eIF4F subunit RNAseq",
         y     = paste0("log2(value)")) +
    stat_compare_means(comparisons = my_comparison1, method = "t.test")
  p1$layers[[2]]$aes_params$textsize <- 5  
  p2 <- plotEIF(EIF.score.tcga.all.tumors) + 
    labs(title = paste0("All tumors n = ", number),
         x     = "eIF4E ratio score",
         y     = paste0("log2(value)")) +
    stat_compare_means(comparisons = my_comparison2, method = "t.test")
  p2$layers[[2]]$aes_params$textsize <- 5
  grid.arrange(p1, p2, ncol=2)
}
plot.EIFandScore.all.tumors()



plot.EIFandScore.each.tissue <- function (m){
  EIF.RNAseq.tcga <- lapply(tcga.study.list, get.EIF.RNAseq.tcga)
  EIF.score.tcga <- lapply(tcga.study.list, get.EIF.score.tcga)
  EIF.RNAseq.GTEx <- EIF.RNAseq.GTEx[, 
                                     c("EIF4A1","EIF4E",
                                       "EIF4G1","EIF4EBP1",
                                       "RPS6KB1","SMTSD")]
  EIF.score.GTEx <- get.EIF.score.GTEx()
  EIF.score.GTEx <- EIF.score.GTEx[, c("EIF4A1score","EIF4Escore",
                                       "EIF4G1score","EIF4EBP1score",
                                       "RPS6KB1score","SMTSD")]
  EIF.RNAseq.GTEx <- EIF.RNAseq.GTEx[EIF.RNAseq.GTEx$SMTSD == m,]
  EIF.score.GTEx <- EIF.score.GTEx[EIF.score.GTEx$SMTSD == m,]
  medianEIF4G1score <- median(EIF.score.GTEx$EIF4G1score)
  medianEIF4EBP1score <- median(EIF.score.GTEx$EIF4EBP1score)
  # tissue$SMTSD <- NULL
  number <- nrow(EIF.RNAseq.GTEx)
  EIF.RNAseq.GTEx.each.tissues <- melt(EIF.RNAseq.GTEx)
  EIF.score.GTEx.each.tissues <- melt(EIF.score.GTEx)
  my_comparison1 <- list( c("EIF4E", "EIF4G1"), 
                          c("EIF4G1", "EIF4EBP1"), 
                          c("EIF4E", "EIF4EBP1"),
                          c("EIF4E", "RPS6KB1"))
  my_comparison2 <- list( c("EIF4Escore", "EIF4G1score"), 
                          c("EIF4G1score", "EIF4EBP1score"), 
                          c("EIF4Escore", "EIF4EBP1score"),
                          c("EIF4Escore", "RPS6KB1score"))
  if (medianEIF4G1score < medianEIF4EBP1score ) {
    par(mfrow=c(1,2))
    p1 <- plotEIF(EIF.RNAseq.GTEx.each.tissues) +     
      labs(title = paste0(m," n = ", number),
           x     = "eIF4F subunit RNAseq",
           y     = paste0("log2(value)")) +
      stat_compare_means(comparisons = my_comparison1, method = "t.test")
    p1$layers[[2]]$aes_params$textsize <- 5  
    p2 <- plotEIF(EIF.score.GTEx.each.tissues) + 
      labs(title = paste0(m," n = ", number),
           x     = "eIF4E ratio score",
           y     = paste0("log2(value)")) +
      stat_compare_means(comparisons = my_comparison2, method = "t.test")
    p2$layers[[2]]$aes_params$textsize <- 5
    grid.arrange(p1, p2, ncol=2)
    print(paste("EIF is inhibited in", m))
  } else {
    print(paste("EIF is activated in", m))
  }
}





























##############################################################
## plot Proteinarray data from all TCGA cancer study groups ##
##############################################################
plot.EIF.tcga <- function(EIF){
  # Get EIF RNAseq data from all TCGA study groups
  tcga.pro.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)", getCancerStudies(mycgds)$name), ]
  # "tcag_study_list" contains all the tcga cancer studies
  tcga.study.list <- tcga.pro.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  # use lappy to pull out all the caselists within tcga.study.list
  # because we named each elements in tcga.study.list,
  # lappy will return a large list, each element (with a cancer study name)
  # in that list is a data-table
  tcga.pro.caselist <- lapply(tcga.study.list, caselist)
  tcga.pro.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  # for example, tcga.pro.caselist[[1]] shows the dataframe of caselist
  # in laml study group.
  # to choose case_list_id that is labeled with laml_tcga_rna_seq_v2_mrna,
  # we use the following tcag_provisional_caselist[[1][8,1]
  # a <- tcga.pro.caselist[[1]][
  # grep("tcga_rna_seq_v2_mrna", tcga.pro.caselist[[1]]$case_list_id),
  # ][1,1]
  # b <- tcga.pro.geneticprofile[[1]][
  # grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
  # tcga.pro.geneticprofile[[1]]$genetic_profile_name), ][1,1]
  # how do we do this for all study groups from [[1]] to  [[32]]?
  caselist.protein <- function(x) {
    tcga.pro.caselist[[x]][
      grep("tcga_rppa",
           tcga.pro.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.protein <- function(x) {
    tcga.pro.geneticprofile[[x]][
      # double backslash \\ suppress the special meaning of ( )
      # in regular expression
      grep("Protein expression Z-scores",
           tcga.pro.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
  # test the functions: caselist.RNAseq () and geneticprofile.RNAseq ()
  # caselist.RNAseq = caselist.RNAseq ('acc_tcga')
  # geneticprofile.RNAseq = geneticprofile.RNAseq ('acc_tcga')
  # Wrap two functions: geneticprofile.RNAseq(x), caselist.RNAseq(x)
  # within TCGA_ProfileData_RNAseq(x)
  tcga.profiledata.protein <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
                   genename,
                   geneticprofile,
                   caselist)
  }
  EIF.tcga.protein <- function(x, y) {
    tcga.profiledata.protein(x,
                             geneticprofile.protein(y),
                             caselist.protein(y))
  }
  EIF.protein.tcga.all <- function(x) {
    test <- lapply(tcga.study.list,
                   function(y) mapply(EIF.tcga.protein, x, y))
    df2 <- melt(test)
    colnames(df2) <- c("RNAseq", "EIFgene", "TCGAstudy")
    df2 <- data.frame(df2)
  }
  df2 <- EIF.protein.tcga.all("SCD")
  df2$EIFgene <- as.factor(df2$EIFgene)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  # plot EIF gene expression across all TCGA groups ##
  m <- paste0("SCD", ".", "SCD")
  mean <- within(df2[df2$EIFgene == m,], # TCGAstudy is one column in df2
                 TCGAstudy <- reorder(TCGAstudy, RNAseq, median))
  print(
    ggplot(mean,
           aes(x     = TCGAstudy,
               y     = log2(RNAseq),
               color = TCGAstudy)) +
      geom_boxplot(alpha    = .01,
                   width    = .5,
                   position = position_dodge(width = .9)) +
      labs(x = "Tumor types (TCGA)",
           y = paste0("log2(", "SCD", " protein level)")) +
      theme(axis.title  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            axis.text.x = element_text(size   = 9,
                                       angle  = 45,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.text.y = element_text(size   = 9,
                                       angle  = 0,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.line.x = element_line(color  = "black"),
            axis.line.y = element_line(color  = "black"),
            panel.grid  = element_blank(),
            strip.text  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            legend.position = "none"))
}

plot.EIF.tcga("SCD")
sapply(EIF.gene, plot.EIF.tcga)

################################################
## Get oncogene mutation data from SKCM group ##
################################################
# note there may be some internal bugs for the data labeled as NaN
# https://github.com/cBioPortal/cgdsr/issues/2
getmutations <- function() {
  mutations <- getProfileData(mycgds,
                              c("BRAF","MYC", "NRAS", "PTEN", "TP53"),
                              "skcm_tcga_mutations",
                              "skcm_tcga_all")
  colnames(mutations) <- paste0(colnames(mutations), '.mutations')
  v <- rownames(mutations)
  # each mutation column contains three types of data:
  # mutation (V600E), NAN (wildtype), NA (not sequenced).
  relabel.mutations <- function(gene) {
    mutations[, gene] <- ifelse(
      mutations[, gene] == "NaN", "Wildtype", "Mutated")
  }
  # use sapply , input as a matrix, and output as a matrix too.
  mutations <- sapply(colnames(mutations), relabel.mutations)
  mutations <- as.data.frame(mutations)
  # the sapply function return a new matrix lacking row names.
  # add row names back with the following two lines
  mutations2 <- cbind(Row.Names = v, mutations)
  # mutations <- mutations2[ , -1]
  rownames(mutations) <- mutations2[ , 1]
  return(mutations)
}

# mutations.data <- getmutations("BRAF")
mutations.data <- getmutations()
# mutations.list <- c("BRAF", "NRAS", "AKT", "TP53")

###########################################
## Get oncogene CNV data from SKCM group ##
###########################################
getCNV <- function(x) {
  CNV <- getProfileData(mycgds,
                        x,
                        "skcm_tcga_gistic",
                        "skcm_tcga_all")
  v <- rownames(CNV)
  colnames(CNV) <- paste0(colnames(CNV), '.CNV')
  # change CNV column from numeric to factor
  factor.CNV <- function(gene) {
    CNV[, gene] <- as.factor(CNV[, gene])
  }
  # use sapply, input as a matrix, and output as a matrix too.
  CNV <- sapply(colnames(CNV), factor.CNV)
  CNV <- as.data.frame(CNV)
  CNV2 <- cbind(Row.Names = v, CNV)
  rownames(CNV) <- CNV2[ , 1]
  #  print(attributes(CNV))
  return(CNV)
}

# CNV.data <- getCNV("BRAF")
CNV.data <- getCNV(c("MYC", "NRAS", "PTEN"))

##########################################################################
## examine the correlation between mutation status and EIF RNASeq level ##
##########################################################################
plot.mutations.RNAseq <- function(mutations, RNAseq) {
  mutations.EIF.RNAseq <- cbind(mutations.data, EIF.RNAseq.data)
  mutations.EIF.RNAseq <- na.omit(mutations.EIF.RNAseq)
  print(
    ggplot(mutations.EIF.RNAseq,
           #use[ ,genemutations] not $genemutations for variable in a function
           aes(x     = mutations.EIF.RNAseq[ ,mutations],
               y     = log2(mutations.EIF.RNAseq[, RNAseq]),
               color = mutations.EIF.RNAseq[ ,mutations])) +
      geom_boxplot(alpha = .01,
                   width = .5) +
      labs(x = mutations,
           y = paste("log2(", RNAseq, "RNA counts)")) +
      theme(axis.title        = element_text(face   = "bold",
                                             size   = 9,
                                             color  = "black"),
            axis.text         = element_text(size   = 9,
                                             hjust  = 1,
                                             face   = "bold",
                                             color  = "black"),
            axis.line.x       = element_line(color  = "black"),
            axis.line.y       = element_line(color  = "black"),
            panel.grid        = element_blank(),
            strip.text        = element_text(face   = "bold",
                                             size   = 9,
                                             colour = "black"),
            legend.position   = "none") +
      # stack the dots along the y-axis and group them along x-axis
      geom_dotplot(binaxis    = "y",
                   binwidth   = .1,
                   stackdir   = "center",
                   fill       = NA))
}

plot.mutations.RNAseq("NRAS.mutations", "SCD")
sapply(c("PIK3CA.mutations", "CDKN2A.mutations",
         "NOTCH1.mutations", "TP53.mutations"),
       function(x) mapply(plot.mutations.RNAseq,
                          x,
                          c("EIF4A1","EIF4E","EIF4G1","EIF4EBP1","RPS6KB1")))
# or
sapply(c("EIF4A1","EIF4E","EIF4G1","EIF4EBP1","RPS6KB1"),
       function(y) mapply(plot.mutations.RNAseq,
                          c("PIK3CA.mutations", "CDKN2A.mutations",
                            "NOTCH1.mutations", "TP53.mutations"),
                          y))

#######################################################################
##  check the data distribution, then choose the statistics methods  ##
#######################################################################
stats <- function(RNAseq, mutations) {
  mutations.EIF.RNAseq <- cbind(mutations.data, EIF.RNAseq.data)
  mutations.EIF.RNAseq <- na.omit(mutations.EIF.RNAseq)
  print(
    ggplot(mutations.EIF.RNAseq,
           aes(x      = mutations.EIF.RNAseq[, RNAseq],
               colour = mutations.EIF.RNAseq[, mutations])) +
      geom_density() +
      scale_x_continuous(name = RNAseq) +
      scale_colour_discrete(name = mutations)) # change the legend title
  Mutated <- function(RNAseq, mutations) {
    mutations.EIF.RNAseq[ ,RNAseq][
      mutations.EIF.RNAseq[ ,mutations] == "Mutated"]
  }
  Wildtype <- function(RNAseq, mutations) {
    mutations.EIF.RNAseq[ ,RNAseq][
      mutations.EIF.RNAseq[ ,mutations] == "Wildtype"]
  }
  list <- list(Mutated(RNAseq, mutations),
               Wildtype(RNAseq, mutations))
  nameit <- function(RNAseq, mutations) {
    name <- c(paste(RNAseq, "with", mutations),
              paste(RNAseq, "without", mutations))}
  names(list) <- nameit(RNAseq, mutations)
  result <- lapply(list, shapiro.test)
  a <- t.test(
    mutations.EIF.RNAseq[ ,RNAseq] ~ mutations.EIF.RNAseq[ ,mutations]
  )
  # overwrite the data.name variable of t.test result
  a$data.name <- paste(RNAseq, 'expression by', mutations, 'status')
  b <- wilcox.test(
    mutations.EIF.RNAseq[ ,RNAseq] ~ mutations.EIF.RNAseq[ ,mutations]
  )
  b$data.name <- paste(RNAseq, 'expression by', mutations, 'status')
  print(result)
  print(a)
  print(b)
}

stats("SCD", "BRAF.mutations")
sapply(c("EIF4A1","EIF4E","EIF4G1","EIF4EBP1","RPS6KB1"),
       function(x)
         mapply(stats,
                x,
                c("PIK3CA.mutations", "CDKN2A.mutations",
                  "NOTCH1.mutations", "TP53.mutations")))
# or
sapply(c("PIK3CA.mutations", "CDKN2A.mutations",
         "NOTCH1.mutations", "TP53.mutations"),
       function(y)
         mapply(stats,
                c("EIF4A1","EIF4E","EIF4G1","EIF4EBP1","RPS6KB1"),
                y))

#############################################################################
## compare the correlation between oncogene CNV and EIF gene RNA-seq data ##
#############################################################################
# make a large function to plot all genes
plot.CNV.RNAseq <- function(geneCNV, RNAseq) {
  CNV.EIF.RNAseq <- cbind(CNV.data, EIF.RNAseq.data)
  # na.omit cannot eleminate NaN here!
  toBeRemoved <- which(CNV.EIF.RNAseq$BRAF.CNV == "NaN")
  CNV.EIF.RNAseq <- CNV.EIF.RNAseq[-toBeRemoved,]
  print(
    ggplot(CNV.EIF.RNAseq,
           aes(x     = CNV.EIF.RNAseq[, geneCNV],
               y     = log2(CNV.EIF.RNAseq[, RNAseq]),
               color = CNV.EIF.RNAseq[, geneCNV])) +
      geom_boxplot(alpha      = .01,
                   width      = .5) +
      geom_dotplot(binaxis    = "y",
                   binwidth   = .1,
                   stackdir   = "center",
                   fill       = NA) +
      scale_x_discrete(limits = c("-2", "-1", "0", "1", "2"), # skip NaN data
                       labels = c("Homdel", "Hetlos", "Diploid", "Gain", "Amp"),
                       drop   = FALSE) +
      labs(x = geneCNV,
           y = paste("log2(", RNAseq, "RNA counts)")) +
      theme(axis.title      = element_text(face   = "bold",
                                           size   = 9,
                                           color  = "black"),
            axis.text       = element_text(size   = 9,
                                           hjust  = 1,
                                           face   = "bold",
                                           color  = "black"),
            axis.line.x     = element_line(color  = "black"),
            axis.line.y     = element_line(color  = "black"),
            panel.grid      = element_blank(),
            strip.text      = element_text(face   = "bold",
                                           size   = 9,
                                           colour = "black"),
            legend.position = "none")
  )
  a <- leveneTest(CNV.EIF.RNAseq[, RNAseq] ~ CNV.EIF.RNAseq[, geneCNV])
  a$data.name <- paste(RNAseq,'expression by', geneCNV, 'status')
  b <- fligner.test(CNV.EIF.RNAseq[, RNAseq] ~ CNV.EIF.RNAseq[, geneCNV])
  b$data.name <- paste(RNAseq,'expression by', geneCNV, 'status')
  print(b)
  print(a)
}
# use all combinations of geneCNV and RNAseq
sapply(c("MYC.CNV","NRAS.CNV","PTEN.CNV"),
       function(x)
         mapply(plot.CNV.RNAseq, x, c("EIF4A1","EIF4E","EIF4G1","EIF4EBP1","RPS6KB1")))
# or
sapply(c("BRAF", "NRAS", "PTEN", "SCD", "FASN"),
       function(y)
         mapply(plot.CNV.RNAseq, c("BRAF.CNV","NRAS.CNV", "PTEN.CNV"), y))
# statistic comparision of SCD expression
# between PTEN heterdeletion and diploid
# (pten loss correlates with scd decrease?)

##################################################################
##  Kaplan-Meier curve with clinic and mutation data from SKCM  ##
##################################################################
plot.km.mut.skcm <- function(ge) {
  mycancerstudy <- getCancerStudies(mycgds)[198, 1]
  mycaselist <- getCaseLists(mycgds, mycancerstudy)[4, 1]
  skcm.clinicaldata <- getClinicalData(mycgds, mycaselist)
  skcm.clinicaldata$rn <- rownames(skcm.clinicaldata)
  mutations.data <- getmutations()
  mutations.data <- na.omit(mutations.data)
  mutations.data$rn <- rownames(mutations.data)
  df <- join_all(list(skcm.clinicaldata[c("OS_MONTHS",
                                          "OS_STATUS",
                                          "rn")],
                      mutations.data),
                 by   = "rn",
                 type = "full")
  df <- as.data.frame(df)
  df <- na.omit(df)
  df$SurvObj <- with(df, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))
  fit <- function(x) {
    survfit(SurvObj ~ df[[x]], data = df, conf.type = "log-log")
  }
  km <- fit(ge)
  black.bold.12pt <- element_text(face   = "bold",
                                  size   = 12,
                                  colour = "black")
  print(
    ggplot2::autoplot(km,
             xlab = "Months",
             ylab = "Survival Probability",
             main = paste("Kaplan-Meier plot", ge)) +
      theme(axis.title      = black.bold.12pt,
            axis.text       = black.bold.12pt,
            axis.line.x     = element_line(color  = "black"),
            axis.line.y     = element_line(color  = "black"),
            panel.grid      = element_blank(),
            strip.text      = black.bold.12pt,
            legend.text     = black.bold.12pt ,
            legend.title    = black.bold.12pt ,
            legend.justification = c(1,1)) +
      scale_fill_discrete(name = ge))
  stats <- function(x) {
    survdiff(SurvObj ~ df[[x]], data = df, rho = 1)
  }
  print(ge)
  print(stats(ge))
}

plot.km.mut.skcm("BRAF.mutations")
mutation.list <- c("BRAF.mutations",
                   "NRAS.mutations",
                   "PTEN.mutations",
                   "TP53.mutations")
names(mutation.list) <- mutation.list
sapply(mutation.list, plot.km.mut.skcm)

####################################################################
##  Kaplan-Meier curve with clinic and EIF RNASeq data from SKCM  ##
####################################################################
plot.km.EIF.skcm <- function(EIF) {
  mycancerstudy <- getCancerStudies(mycgds)[217, 1]        # "laml_tcga"
  mycaselist <- getCaseLists(mycgds, mycancerstudy)[4, 1]  # "hnsc_tcga_all"
  skcm.clinicaldata <- getClinicalData(mycgds, mycaselist)
  skcm.clinicaldata$rn <- rownames(skcm.clinicaldata)
  skcm.RNAseq.data <- getProfileData(mycgds,
                                     EIF,
                                     "tgct_tcga_rna_seq_v2_mrna",
                                     "tgct_tcga_all")
  skcm.RNAseq.data <- as.data.frame(skcm.RNAseq.data)
  skcm.RNAseq.data$rn <- rownames(skcm.RNAseq.data)
  df <- join_all(list(skcm.clinicaldata[c("OS_MONTHS", "OS_STATUS", "rn")],
                      skcm.RNAseq.data),
                 by   = "rn",
                 type = "full")
  df <- na.omit(df)
  number <- nrow(df)
  df$Group[df[[EIF]] < quantile(df[[EIF]], prob = 0.2)] = "Bottom 20%"
  df$Group[df[[EIF]] > quantile(df[[EIF]], prob = 0.8)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))
  df <- na.omit(df)
  km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 0)
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
  black.bold.12pt <- element_text(face   = "bold",
                                  size   = 12,
                                  colour = "black")
  print(
    ggplot2::autoplot(km,
             xlab = "Months",
             ylab = "Survival Probability",
             main = paste0("Kaplan-Meier plot of TCGA Sarcoma group (", 
                           number," cases)")) +
      theme(axis.title           = black.bold.12pt,
            axis.text            = black.bold.12pt,
            axis.line.x          = element_line(color  = "black"),
            axis.line.y          = element_line(color  = "black"),
            panel.grid           = element_blank(),
            strip.text           = black.bold.12pt,
            legend.text          = black.bold.12pt ,
            legend.title         = black.bold.12pt ,
            legend.position      = c(1,1),
            legend.justification = c(1,1)) +
      guides(fill = FALSE) +
      scale_color_manual(values = c("red", "blue"),
                         name   = paste(EIF, "mRNA expression"),
                         breaks = c("Bottom 20%", "Top 20%"),
                         labels = c("Bottom 20%, n = 53",
                                    "Top 20%, n = 53")) +
      geom_point(size = 0.25) +
      annotate("text",
               x     = 200,
               y     = 0.8,
               label = paste("log-rank test, p.val = ", p.val),
               size  = 4.5,
               hjust = 1,
               fontface = "bold"))
  # rho = 1 the Gehan-Wilcoxon test
  print(EIF)
  print(stats)
  #  fit = survfit(SurvObj ~ df$Group, data = df)
  #  tst <- comp(fit)$tests$lrTests
  #  print(tst)
}

plot.km.EIF.skcm("EIF4G1")
sapply(EIF.gene, plot.km.EIF.skcm)

#########################################################################################
## Kaplan-Meier curve with clinic and EIF RNAseq data from all TCGA provisional groups ##
#########################################################################################
plot.km.all.pro.tcga <- function(EIF) {
  #  mycgds <- CGDS("http://www.cbioportal.org/")
  #  test(mycgds)
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  tcga.pro.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)",
         getCancerStudies(mycgds)$name), ]
  tcga.study.list <- tcga.pro.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  tcga.pro.caselist <- lapply(tcga.study.list,
                              caselist)
  tcga.pro.geneticprofile <- lapply(tcga.study.list,
                                    geneticprofile)
  caselist.RNAseq <- function(x) {
    tcga.pro.caselist[[x]][grep("tcga_rna_seq_v2_mrna",
                                tcga.pro.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.RNAseq <- function(x) {
    tcga.pro.geneticprofile[[x]][
      grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
           tcga.pro.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
  tcga.profiledata.RNAseq <- function(genename,
                                      geneticprofile,
                                      caselist) {
    getProfileData(mycgds,
                   genename,
                   geneticprofile,
                   caselist)
  }
  tcga.gene.RNAseq <- function(x, y) {
    tcga.profiledata.RNAseq(x,
                            geneticprofile.RNAseq(y),
                            caselist.RNAseq(y))
  }
  tcga.EIF.RNAseq <- function(y) {
    tcga.gene.RNAseq(x = EIF, y)
  }
  ## test1 <- SCD.tcga.RNAseq(y = "skcm_tcga")
  ## test ## try to keep patient ID in the rowname
  test2 <- lapply(tcga.study.list, tcga.EIF.RNAseq)
  for(x in 1:32)
  {
    test2[[x]]$case.id <- rownames(test2[[x]])
    message("test2 = ", x)
  }
  df1 <- melt(test2)
  colnames(df1) <- c("case.id",
                     "EIFgene",
                     "RNAseq",
                     "TCGAstudy")
  all.tcga.EIF.RNAseq <- data.frame(df1)
  all.tcga.EIF.RNAseq$TCGAstudy <- as.factor(all.tcga.EIF.RNAseq$TCGAstudy)
  message("RNAseq data retrieved")
  ##### retrieve clinic data from all tcga groups #####
  tcga.clinic.data <- function(x) {
    print(x)
    url <- function(x){
      url <- "http://www.cbioportal.org/webservice.do?cmd=getClinicalData&case_set_id="
      url <- paste0(url, x, "_all")
      return(url)
    }
    # testurl <- url("acc_tcga")
    # tesereq <- GET(url("acc_tcga"))
    req <- function(x) {GET(url(x))}
    # req <- req("acc_tcga")
    clinical_data <- function(x) {content(req(x),
                                          type      = 'text/tab-separated-values',
                                          col_names = T,
                                          col_types = NULL)}
    data <- clinical_data(x)
    data <- data[c("OS_MONTHS",
                   "OS_STATUS",
                   "CASE_ID")]
  }
  # three datasets donot have OS data and cause bugs remove them
  bug.data.set <- names(tcga.study.list) %in% c("meso_tcga", "pcpg_tcga", "ucs_tcga")
  tcga.study.list <- tcga.study.list[!bug.data.set]
  all.tcga.clinic.data <- lapply(tcga.study.list, tcga.clinic.data)
  all.tcga.clinic.data <- melt(all.tcga.clinic.data)
  all.tcga.clinic.data <- all.tcga.clinic.data[c("OS_STATUS",
                                                 "CASE_ID",
                                                 "value",
                                                 "L1")]
  colnames(all.tcga.clinic.data) <- c("OS_STATUS",
                                      "case.id",
                                      "OS_MONTHS",
                                      "TCGAstudy")
  all.tcga.clinic.data$case.id <- str_replace_all(all.tcga.clinic.data$case.id,
                                                  '-',
                                                  '.')
  message("clinical data retrieved")
  df <- join_all(list(all.tcga.clinic.data[c("OS_MONTHS",
                                             "OS_STATUS",
                                             "case.id")],
                      all.tcga.EIF.RNAseq[c("case.id",
                                            "RNAseq",
                                            "TCGAstudy")]),
                 by   = "case.id",
                 type = "full")
  df <- na.omit(df)
  message("clinical and RNAseq data combined")
  df$Group[df$RNAseq < quantile(df$RNAseq, prob = 0.2)] = "Bottom 20%"
  df$Group[df$RNAseq > quantile(df$RNAseq, prob = 0.8)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))
  #  df <- na.omit(df)
  km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 0)
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
  black.bold.12pt <- element_text(face   = "bold",
                                  size   = 12,
                                  colour = "black")
  print(
    ggplot2::autoplot(km,
             xlab = "Months",
             ylab = "Survival Probability",
             main = paste("Kaplan-Meier plot", EIF, 
                          "RNA expression in all TCGA provisional groups")) +
      theme(axis.title           = black.bold.12pt,
            axis.text            = black.bold.12pt,
            axis.line.x          = element_line(color  = "black"),
            axis.line.y          = element_line(color  = "black"),
            panel.grid           = element_blank(),
            strip.text           = black.bold.12pt,
            legend.text          = black.bold.12pt ,
            legend.title         = black.bold.12pt ,
            legend.position      = c(1,1),
            legend.justification = c(1,1)) +
      guides(fill = FALSE) +
      scale_color_manual(values = c("red", "blue"),
                         name   = paste(EIF, "mRNA expression"),
                         breaks = c("Bottom 20%", "Top 20%"),
                         labels = c("Bottom 20%, n = 1859",
                                    "Top 20%, n = 1859")) +
      geom_point(size = 0.25) +
      annotate("text",
               x     = 300,
               y     = 0.85,
               label = paste("log-rank test, p.val = ", p.val),
               size  = 4.5,
               hjust = 1,
               fontface = "bold")
  )
  # rho = 1 the Gehan-Wilcoxon test
  #  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 1)
  print(EIF)
  print(stats)
}

plot.km.all.pro.tcga("EIF4G1")
sapply(EIF.gene, plot.km.all.pro.tcga)

#######################################################################################
## Kaplan-Meier curve with clinic and EIF RNAseq data from all TCGA pancancer groups ##
#######################################################################################
plot.km.all.pan.tcga <- function(EIF) {
  #  mycgds <- CGDS("http://www.cbioportal.org/")
  #  test(mycgds)
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  pan.tcga.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)",
         getCancerStudies(mycgds)$name), ]
  pan.tcga.study.list <- pan.tcga.studies$cancer_study_id
  names(pan.tcga.study.list) <- pan.tcga.study.list
  pan.tcga.caselist <- lapply(pan.tcga.study.list,
                              caselist)
  pan.tcga.geneticprofile <- lapply(pan.tcga.study.list,
                                    geneticprofile)
  pan.caselist.RNAseq <- function(x) {
    pan.tcga.caselist[[x]][grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
                                pan.tcga.caselist[[x]]$case_list_id), ][1, 1]
  } # pancancer group does not contain OS data
  pan.geneticprofile.RNAseq <- function(x) {
    pan.tcga.geneticprofile[[x]][
      grep("mRNA Expression, RSEM",
           pan.tcga.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
  tcga.profiledata.RNAseq <- function(genename,
                                      geneticprofile,
                                      caselist) {
    getProfileData(mycgds,
                   genename,
                   geneticprofile,
                   caselist)
  }
  pan.tcga.gene.RNAseq <- function(x, y) {
    tcga.profiledata.RNAseq(x,
                            pan.geneticprofile.RNAseq(y),
                            pan.caselist.RNAseq(y))
  }
  pan.tcga.EIF.RNAseq <- function(y) {
    pan.tcga.gene.RNAseq(x = EIF, y)
  }
  ## test1 <- SCD.tcga.RNAseq(y = "skcm_tcga")
  ## test ## try to keep patient ID in the rowname
  test2 <- lapply(pan.tcga.study.list, pan.tcga.EIF.RNAseq)
  for(x in 1:32)
  {
    test2[[x]]$case.id <- rownames(test2[[x]])
    message("test2 = ", x)
  }
  df1 <- melt(test2)
  colnames(df1) <- c("case.id",
                     "EIFgene",
                     "RNAseq",
                     "TCGAstudy")
  all.pan.tcga.EIF.RNAseq <- data.frame(df1)
  all.pan.tcga.EIF.RNAseq$TCGAstudy <- as.factor(all.pan.tcga.EIF.RNAseq$TCGAstudy)
  message("RNAseq data retrieved")
  ##### retrieve clinic data from all tcga groups #####
  tcga.clinic.data <- function(x) {
    print(x)
    url <- function(x){
      url <- "http://www.cbioportal.org/webservice.do?cmd=getClinicalData&case_set_id="
      url <- paste0(url, x, "_all")
      return(url)
    }
    # testurl <- url("acc_tcga")
    # tesereq <- GET(url("acc_tcga"))
    req <- function(x) {GET(url(x))}
    # req <- req("acc_tcga")
    clinical_data <- function(x) {content(req(x),
                                          type      = 'text/tab-separated-values',
                                          col_names = T,
                                          col_types = NULL)}
    data <- clinical_data(x)
    data <- data[c("OS_MONTHS",
                   "OS_STATUS",
                   "CASE_ID")]
  }
  
  pro.tcga.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)", getCancerStudies(mycgds)$name), ]
  # "pro.tcga.study.list" contains all the tcga provisional cancer studies
  pro.tcga.study.list <- pro.tcga.studies$cancer_study_id
  names(pro.tcga.study.list) <- pro.tcga.study.list
  # three datasets donot have OS data and cause bugs remove them
  bug.data.set <- names(pro.tcga.study.list) %in% c("meso_tcga", "pcpg_tcga", "ucs_tcga")
  pro.tcga.study.list <- pro.tcga.study.list[!bug.data.set]
  all.tcga.clinic.data <- lapply(pro.tcga.study.list, tcga.clinic.data)
  all.tcga.clinic.data <- melt(all.tcga.clinic.data)
  all.tcga.clinic.data <- all.tcga.clinic.data[c("OS_STATUS",
                                                 "CASE_ID",
                                                 "value",
                                                 "L1")]
  colnames(all.tcga.clinic.data) <- c("OS_STATUS",
                                      "case.id",
                                      "OS_MONTHS",
                                      "TCGAstudy")
  all.tcga.clinic.data$case.id <- str_replace_all(all.tcga.clinic.data$case.id,
                                                  '-',
                                                  '.')
  message("clinical data retrieved")
  df <- join_all(list(all.tcga.clinic.data[c("OS_MONTHS",
                                             "OS_STATUS",
                                             "case.id")],
                      all.pan.tcga.EIF.RNAseq[c("case.id",
                                                "RNAseq",
                                                "TCGAstudy")]),
                 by   = "case.id",
                 type = "full")
  df <- na.omit(df)
  message("clinical and RNAseq data combined")
  number <- nrow(df)
  df$Group[df$RNAseq < quantile(df$RNAseq, prob = 0.2)] = "Bottom 20%"
  df$Group[df$RNAseq > quantile(df$RNAseq, prob = 0.8)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))
  #  df <- na.omit(df)
  km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 0)
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
  black.bold.12pt <- element_text(face   = "bold",
                                  size   = 12,
                                  colour = "black")
  p <- ggplot2::autoplot(km,
                xlab = "Months",
                ylab = "Survival Probability",
                main = paste0("Kaplan-Meier plot in TCGA pancancer groups (", number," cases)"),
                xlim = c(0, 300)) +
    #      scale_fill_discrete(name="Experimental\nCondition",
    #                          breaks=c("Top 20%", "Bottom 20%"),
    #                          labels=c("Top 20%, n = 1931", "Bottom 20%, n = 1931"))+
    theme(axis.title           = black.bold.12pt,
          axis.text            = black.bold.12pt,
          axis.line.x          = element_line(color  = "black"),
          axis.line.y          = element_line(color  = "black"),
          panel.grid           = element_blank(),
          strip.text           = black.bold.12pt,
          legend.text          = black.bold.12pt,
          legend.title         = black.bold.12pt,
          legend.position      = c(1,1), #place the top-right corner of the legend
          legend.justification = c(1,1)) + #in the top-right corner fot eh graph
    guides(fill = FALSE) +
    scale_color_manual(values = c("red", "blue"),
                       name   = paste(EIF, "mRNA expression"),
                       breaks = c("Bottom 20%", "Top 20%"),
                       labels = c("Bottom 20%, n = 1931",
                                  "Top 20%, n = 1931")) +
    geom_point(size = 0.25) +
    annotate("text",
             x     = 300,
             y     = 0.90,
             label = paste("log-rank test, p.val = ", p.val),
             size  = 4.5,
             hjust = 1,
             fontface = "bold")
  
  # rho = 1 the Gehan-Wilcoxon test
  #  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 1)
  #  my_text <- "This text is at x=0.7 and y=0.8!"
  #  my_grob = grid.text(my_text, x=0.7,  y=0.8, gp=gpar(col="firebrick", fontsize=14, fontface="bold"))
  #  p + annotation_custom(my_grob)
  print(p)
  print(stats)
}

plot.km.all.pan.tcga("EIF4E")
sapply(EIF.gene, plot.km.all.pan.tcga)
