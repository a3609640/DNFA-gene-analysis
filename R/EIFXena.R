## the following script perform PCA on RNA-Seq data of seven DNFA genes 
## from SKCM amd GTEX dataset with R package "ggfortify".
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(gridExtra)
library(reshape2)
library(survival)
library(survMisc)
library(survminer)

## read.csv will transform characters into factors  
get.EIF.TCGA.GTEX.RNAseq.long <- function () {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data", "EIFTCGAGTEX2.csv"), 
                            header = TRUE, sep = ",")
  EIF.TCGA.GTEX.RNAseq.long <- melt(EIF.TCGA.GTEX[, 1:17])
  colnames(EIF.TCGA.GTEX.RNAseq.long) <- c("sample", "study", 
                                           "sample.type", "primary.disease", 
                                           "variable", "value")
  EIF.TCGA.GTEX.RNAseq.long <- EIF.TCGA.GTEX.RNAseq.long[
    EIF.TCGA.GTEX.RNAseq.long$value != 0,]
  EIF.TCGA.GTEX.RNAseq.long <- na.omit(EIF.TCGA.GTEX.RNAseq.long)
  tumor.type <- c("Metastatic", "Primary Tumor", 
                  "Recurrent Tumor", "Solid Tissue Normal", 
                  "Normal Tissue", "Cell Line")
  EIF.TCGA.GTEX.RNAseq.long <- EIF.TCGA.GTEX.RNAseq.long[
    EIF.TCGA.GTEX.RNAseq.long$sample.type %in% tumor.type,]
  EIF.TCGA.GTEX.RNAseq.long <- droplevels(EIF.TCGA.GTEX.RNAseq.long)
  EIF.TCGA.GTEX.RNAseq.long$sample.type <- factor(
    EIF.TCGA.GTEX.RNAseq.long$sample.type, levels = tumor.type)
  return(EIF.TCGA.GTEX.RNAseq.long)
}

##
get.DNFA.TCGA.GTEX.RNAseq.long <- function () {
  DNFA.TCGA.GTEX <- read.csv(file.path("project-data", "DNFASKCMandGTEX.csv"), 
                            header = TRUE, sep = ",")
  DNFA.TCGA.GTEX.RNAseq.long <- melt(DNFA.TCGA.GTEX[, 1:14])
  colnames(DNFA.TCGA.GTEX.RNAseq.long) <- c("sample", "study", 
                                            "sample.type","primary.disease", 
                                            "variable", "value")
  DNFA.TCGA.GTEX.RNAseq.long <- DNFA.TCGA.GTEX.RNAseq.long[
    DNFA.TCGA.GTEX.RNAseq.long$value != 0,]
  DNFA.TCGA.GTEX.RNAseq.long <- na.omit(DNFA.TCGA.GTEX.RNAseq.long)
  tumor.type <- c("Metastatic", "Primary Tumor", 
                  "Recurrent Tumor", "Solid Tissue Normal", 
                  "Normal Tissue", "Cell Line")
  DNFA.TCGA.GTEX.RNAseq.long <- DNFA.TCGA.GTEX.RNAseq.long[
    DNFA.TCGA.GTEX.RNAseq.long$sample.type %in% tumor.type,]
  DNFA.TCGA.GTEX.RNAseq.long <- droplevels(DNFA.TCGA.GTEX.RNAseq.long)
  DNFA.TCGA.GTEX.RNAseq.long$sample.type <- factor(DNFA.TCGA.GTEX.RNAseq.long$sample.type, 
                                                  levels = tumor.type)
  return(DNFA.TCGA.GTEX.RNAseq.long)
}

##
get.EIF.TCGA.RNAseq.long <- function () {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data", 
                                      "EIFTCGAGTEX.csv"), 
                            header = TRUE, sep = ",")
  EIF.TCGA.GTEX.RNAseq.long <- melt(EIF.TCGA.GTEX[, 1:10])
  colnames(EIF.TCGA.GTEX.RNAseq.long) <- c("sample", "study", 
                                           "sample.type", "primary.disease", 
                                           "variable", "value")
  EIF.TCGA.RNAseq.long <- EIF.TCGA.GTEX.RNAseq.long[
    EIF.TCGA.GTEX.RNAseq.long$study == 'TCGA',]
  EIF.TCGA.RNAseq.long <- EIF.TCGA.RNAseq.long[
    EIF.TCGA.RNAseq.long$value != 0,]
  tumor.type <- c("Metastatic","Primary Tumor",
                  "Recurrent Tumor","Normal Tissue",
                  "Solid Tissue Normal")
  EIF.TCGA.RNAseq.long <- EIF.TCGA.RNAseq.long[
    EIF.TCGA.RNAseq.long$sample.type %in% tumor.type,]
  EIF.TCGA.RNAseq.long <- droplevels(EIF.TCGA.RNAseq.long)
  return(EIF.TCGA.RNAseq.long)
}

##
get.EIF.GTEX.RNAseq.long <- function () {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data", 
                                      "EIFTCGAGTEX.csv"), 
                            header = TRUE, sep = ",")
  EIF.TCGA.GTEX.RNAseq.long <- melt(EIF.TCGA.GTEX[, 1:10])
  colnames(EIF.TCGA.GTEX.RNAseq.long) <- c("sample", "study", 
                                           "sample.type", "primary.disease", 
                                           "variable", "value")
  EIF.GTEX.RNAseq.long <- EIF.TCGA.GTEX.RNAseq.long[
    EIF.TCGA.GTEX.RNAseq.long$study == 'GTEX',]
  EIF.GTEX.RNAseq.long <- EIF.GTEX.RNAseq.long[
    EIF.GTEX.RNAseq.long$value != 0,]
  EIF.GTEX.RNAseq.long <- na.omit(EIF.GTEX.RNAseq.long)
#  EIF.GTEX.RNAseq.long <- EIF.GTEX.RNAseq.long[EIF.GTEX.RNAseq.long$sample.type %in% tumor.type,]
  EIF.GTEX.RNAseq.long <- droplevels(EIF.GTEX.RNAseq.long)
  return(EIF.GTEX.RNAseq.long)
}

##
get.EIF.TCGA.GTEX.score.long <- function () {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data", 
                                      "EIFTCGAGTEX.csv"), 
                            header = TRUE, sep = ",")
  EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX
  EIF.TCGA.GTEX.score$EIF4A1 <- (EIF.TCGA.GTEX$EIF4A1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4A1")] <- "EIF4A1:EIF4E ratio" 
  EIF.TCGA.GTEX.score$EIF4G1 <- (EIF.TCGA.GTEX$EIF4G1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4G1")] <- "EIF4G1:EIF4E ratio" 
  EIF.TCGA.GTEX.score$EIF4EBP1 <- (EIF.TCGA.GTEX$EIF4EBP1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4EBP1")] <- "EIF4EBP1:EIF4E ratio" 
  EIF.TCGA.GTEX.score$RPS6KB1 <- (EIF.TCGA.GTEX$RPS6KB1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "RPS6KB1")] <- "RPS6KB1:EIF4E ratio" 
  EIF.TCGA.GTEX.score$MYC <- (EIF.TCGA.GTEX$MYC - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "MYC")] <- "MYC:EIF4E ratio" 
  EIF.TCGA.GTEX.score$EIF4E <- (EIF.TCGA.GTEX$EIF4E - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4E")] <- "EIF4E:EIF4E ratio" 
  EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.GTEX.score[, 1:10])
  colnames(EIF.TCGA.GTEX.score.long) <- c("sample", "study", "sample.type", 
                                          "primary.disease", "variable", "value")
  tumor.type <- c("Metastatic", "Primary Tumor", "Recurrent Tumor",
                  "Solid Tissue Normal", "Normal Tissue", "Cell Line")
  EIF.TCGA.GTEX.score.long <- EIF.TCGA.GTEX.score.long[
    EIF.TCGA.GTEX.score.long$sample.type %in% tumor.type,]
  EIF.TCGA.GTEX.score.long <- droplevels(EIF.TCGA.GTEX.score.long)
  EIF.TCGA.GTEX.score.long$sample.type <- factor(EIF.TCGA.GTEX.score.long$sample.type, 
                                                 levels = tumor.type)
  return(EIF.TCGA.GTEX.score.long)
}

##
get.EIF.TCGA.score.long <- function () {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data", 
                                      "EIFTCGAGTEX.csv"), 
                            header = TRUE, sep = ",")
  EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX
  EIF.TCGA.GTEX.score$EIF4A1 <- (EIF.TCGA.GTEX$EIF4A1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4A1")] <- "EIF4A1:EIF4E ratio" 
  EIF.TCGA.GTEX.score$EIF4G1 <- (EIF.TCGA.GTEX$EIF4G1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4G1")] <- "EIF4G1:EIF4E ratio" 
  EIF.TCGA.GTEX.score$EIF4EBP1 <- (EIF.TCGA.GTEX$EIF4EBP1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4EBP1")] <- "EIF4EBP1:EIF4E ratio" 
  EIF.TCGA.GTEX.score$RPS6KB1 <- (EIF.TCGA.GTEX$RPS6KB1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "RPS6KB1")] <- "RPS6KB1:EIF4E ratio" 
  EIF.TCGA.GTEX.score$MYC <- (EIF.TCGA.GTEX$MYC - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "MYC")] <- "MYC:EIF4E ratio" 
  EIF.TCGA.GTEX.score$EIF4E <- (EIF.TCGA.GTEX$EIF4E - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4E")] <- "EIF4E:EIF4E ratio" 
  EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.GTEX.score[, 1:10])
  colnames(EIF.TCGA.GTEX.score.long) <- c("sample", "study", "sample.type", 
                                           "primary.disease", "variable", "value")
  EIF.TCGA.score.long <- EIF.TCGA.GTEX.score.long[EIF.TCGA.GTEX.score.long$study == 'TCGA',]
  tumor.type <- c("Metastatic","Primary Tumor",
                  "Recurrent Tumor","Normal Tissue",
                  "Solid Tissue Normal")
  EIF.TCGA.score.long <- EIF.TCGA.score.long[EIF.TCGA.score.long$sample.type %in% tumor.type,]
  EIF.TCGA.score.long <- droplevels(EIF.TCGA.score.long)
  return(EIF.TCGA.score.long)
}

##
get.EIF.GTEX.score.long <- function () {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data", 
                                      "EIFTCGAGTEX.csv"), 
                            header = TRUE, sep = ",")
  EIF.TCGA.GTEX.score <- EIF.TCGA.GTEX
  EIF.TCGA.GTEX.score$EIF4A1 <- (EIF.TCGA.GTEX$EIF4A1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4A1")] <- "EIF4A1:EIF4E ratio" 
  EIF.TCGA.GTEX.score$EIF4G1 <- (EIF.TCGA.GTEX$EIF4G1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4G1")] <- "EIF4G1:EIF4E ratio" 
  EIF.TCGA.GTEX.score$EIF4EBP1 <- (EIF.TCGA.GTEX$EIF4EBP1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4EBP1")] <- "EIF4EBP1:EIF4E ratio" 
  EIF.TCGA.GTEX.score$RPS6KB1 <- (EIF.TCGA.GTEX$RPS6KB1 - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "RPS6KB1")] <- "RPS6KB1:EIF4E ratio" 
  EIF.TCGA.GTEX.score$MYC <- (EIF.TCGA.GTEX$MYC - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "MYC")] <- "MYC:EIF4E ratio" 
  EIF.TCGA.GTEX.score$EIF4E <- (EIF.TCGA.GTEX$EIF4E - EIF.TCGA.GTEX$EIF4E)
  colnames(EIF.TCGA.GTEX.score)[which(names(EIF.TCGA.GTEX.score) == "EIF4E")] <- "EIF4E:EIF4E ratio" 
  EIF.TCGA.GTEX.score.long <- melt(EIF.TCGA.GTEX.score[, 1:10])
  colnames(EIF.TCGA.GTEX.score.long) <- c("sample", "study", "sample.type", 
                                          "primary.disease", "variable", "value")
  EIF.GTEX.score.long <- EIF.TCGA.GTEX.score.long[EIF.TCGA.GTEX.score.long$study == 'GTEX',]
  EIF.GTEX.score.long <- na.omit(EIF.GTEX.score.long)
  EIF.GTEX.score.long <- droplevels(EIF.GTEX.score.long)
  return(EIF.GTEX.score.long)
}

##
plotEIF.TCGA.GTEX <-  function (x) {
  name <- deparse(substitute(x))
  metastatic.number <- nrow(x[x$sample.type == "Metastatic",])
  primary.tumor.number <- nrow(x[x$sample.type == "Primary Tumor",])
  recurrent.tumor.number <- nrow(x[x$sample.type == "Recurrent Tumor",])
  solid.tissue.normal.number <- nrow(x[x$sample.type == "Solid Tissue Normal",])
  cell.line.number <- nrow(x[x$sample.type == "Cell Line",])
  normal.tissue.number <- nrow(x[x$sample.type == "Normal Tissue",])
  black_bold_tahoma_12 <- element_text(color  = "black", 
                                       face   = "bold",
                                       family = "Tahoma", 
                                       size   = 9)
  black_bold_tahoma_12_45 <- element_text(color  = "black",
                                          face   = "bold",
                                          family = "Tahoma", 
                                          size   = 9, 
                                          angle  = 45,
                                          hjust  = 1)
  p1 <- ggplot(data = x,
         aes(x     = sample.type, 
             y     = value, 
             color = sample.type)) +
    facet_grid( ~ variable, 
                scales = "free", 
                space  = "free")+   
    facet_wrap( ~ variable, ncol = 6)+
    geom_violin(trim = FALSE) +
    geom_boxplot(alpha    = .01, 
                 size     = .75,
                 width    = .5,
                 position = position_dodge(width = .9)) +
    labs(x = "sample type",
         y = paste("log2(TPM)")) +
    scale_x_discrete(labels = c("Metastatic"          = paste("Metastatic \n n= ",
                                                              metastatic.number), 
                                "Primary Tumor"       = paste("Primary Tumor \n n= ", 
                                                              primary.tumor.number),
                                "Recurrent Tumor"     = paste("Recurrent Tumor \n n= ", 
                                                              recurrent.tumor.number),
                                "Solid Tissue Normal" = paste("Solid Tissue Normal \n n= ", 
                                                              solid.tissue.normal.number),
                                "Normal Tissue"       = paste("Normal Tissue \n n= ", 
                                                              normal.tissue.number),
                                "Cell Line"           = paste("Cell Line \n n= ", 
                                                              cell.line.number))) +
    theme_bw() +
    theme(plot.title      = black_bold_tahoma_12,
          axis.title      = black_bold_tahoma_12,
          axis.text.x     = black_bold_tahoma_12_45,
          axis.text.y     = black_bold_tahoma_12,
          axis.line.x     = element_line(color = "black"),
          axis.line.y     = element_line(color = "black"),
          panel.grid      = element_blank(),
          legend.position = "none",
          strip.text      = black_bold_tahoma_12) 
  p1 <- p1 + stat_compare_means(method = "anova")
  p2 <- ggplot(data = x,
               aes(x     = variable, 
                   y     = value, 
                   color = variable)) +
    facet_grid( ~ sample.type, 
                scales = "free", 
                space  = "free")+   
    facet_wrap( ~ sample.type, ncol = 6)+
    geom_violin(trim = FALSE) +
    geom_boxplot(alpha    = .01, 
                 size     = .75,
                 width    = .5,
                 position = position_dodge(width = .9)) +
    labs(x = "sample type",
         y = paste("log2(value)")) +
    theme_bw() +
    theme(plot.title      = black_bold_tahoma_12,
          axis.title      = black_bold_tahoma_12,
          axis.text.x     = black_bold_tahoma_12_45,
          axis.text.y     = black_bold_tahoma_12,
          axis.line.x     = element_line(color = "black"),
          axis.line.y     = element_line(color = "black"),
          panel.grid      = element_blank(),
          legend.position = "none",
          strip.text      = black_bold_tahoma_12)
  p2 <- p2 + stat_compare_means(method = "anova")
  print(p1)
  print(p2)
}

##
plotEIF.RNAseq.TCGA <-  function (x) {
  name <- deparse(substitute(x))
  metastatic.number <- nrow(x[x$sample.type == "Metastatic",])
  primary.tumor.number <- nrow(x[x$sample.type == "Primary Tumor",])
  recurrent.tumor.number <- nrow(x[x$sample.type == "Recurrent Tumor",])
  solid.tissue.normal.number <- nrow(x[x$sample.type == "Solid Tissue Normal",])
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
  p1 <- ggplot(data = x,
               aes(x     = sample.type, 
                   y     = value, 
                   color = sample.type)) +
    facet_grid( ~ variable, 
                scales = "free", 
                space  = "free")+   
    facet_wrap( ~ variable, ncol = 6)+
    geom_violin(trim = FALSE) +
    geom_boxplot(alpha    = .01, 
                 size     = .75,
                 width    = .5,
                 position = position_dodge(width = .9)) +
    labs(x = "sample type",
         y = paste("log2(RNA counts)")) +
    scale_x_discrete(labels = c("Metastatic"          = paste("Metastatic \n n= ",
                                                              metastatic.number), 
                                "Primary Tumor"       = paste("Primary Tumor \n n= ", 
                                                              primary.tumor.number),
                                "Recurrent Tumor"     = paste("Recurrent Tumor \n n= ", 
                                                              recurrent.tumor.number),
                                "Solid Tissue Normal" = paste("Solid Tissue Normal \n n= ", 
                                                              solid.tissue.normal.number))) +
    theme_bw() +
    theme(plot.title      = black_bold_tahoma_12,
          axis.title      = black_bold_tahoma_12,
          axis.text.x     = black_bold_tahoma_12_45,
          axis.text.y     = black_bold_tahoma_12,
          axis.line.x     = element_line(color = "black"),
          axis.line.y     = element_line(color = "black"),
          panel.grid      = element_blank(),
          legend.position = "none",
          strip.text      = black_bold_tahoma_12) +
    stat_compare_means(comparisons = list(c("Metastatic", 
                                            "Solid Tissue Normal"), 
                                          c("Primary Tumor", 
                                            "Solid Tissue Normal"), 
                                          c("Recurrent Tumor", 
                                            "Solid Tissue Normal"),
                                          c("Metastatic", 
                                            "Primary Tumor"),
                                          c("Recurrent Tumor", 
                                            "Primary Tumor")), 
                       method = "t.test")  
  
  p2 <- ggplot(data = x,
               aes(x     = variable, 
                   y     = value, 
                   color = variable)) +
    facet_grid( ~ sample.type, 
                scales = "free", 
                space  = "free")+   
    facet_wrap( ~ sample.type, ncol = 6)+
    geom_violin(trim = FALSE) +
    geom_boxplot(alpha    = .01, 
                 size     = .75,
                 width    = .5,
                 position = position_dodge(width = .9)) +
    labs(x = "EIF complex",
         y = paste("log2(RNA counts)")) +
    theme_bw() +
    theme(plot.title      = black_bold_tahoma_12,
          axis.title      = black_bold_tahoma_12,
          axis.text.x     = black_bold_tahoma_12_45,
          axis.text.y     = black_bold_tahoma_12,
          axis.line.x     = element_line(color = "black"),
          axis.line.y     = element_line(color = "black"),
          panel.grid      = element_blank(),
          legend.position = "none",
          strip.text      = black_bold_tahoma_12) +
    stat_compare_means(comparisons = list(c("EIF4A1", "EIF4E"), 
                                          c("EIF4G1", "EIF4E"), 
                                          c("EIF4EBP1", "EIF4E"),
                                          c("RPS6KB1", "EIF4E"),
                                          c("MYC", "EIF4E")), 
                       method = "t.test")  
  print(p1)
  print(p2)
}

##
plotEIF.score.TCGA <-  function (x) {
  name <- deparse(substitute(x))
  metastatic.number <- nrow(x[x$sample.type == "Metastatic",])
  primary.tumor.number <- nrow(x[x$sample.type == "Primary Tumor",])
  recurrent.tumor.number <- nrow(x[x$sample.type == "Recurrent Tumor",])
  solid.tissue.normal.number <- nrow(x[x$sample.type == "Solid Tissue Normal",])
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
  p1 <- ggplot(data = x,
               aes(x     = sample.type, 
                   y     = value, 
                   color = sample.type)) +
    facet_grid( ~ variable, 
                scales = "free", 
                space  = "free")+   
    facet_wrap( ~ variable, 
                ncol = 6)+
    geom_violin(trim = FALSE) +
    geom_boxplot(alpha    = .01, 
                 size     = .75,
                 width    = .5,
                 position = position_dodge(width = .9)) +
    labs(x = "sample type",
         y = paste("log2(value)")) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_x_discrete(labels = c(
      "Metastatic"          = paste("Metastatic \n n= ", 
                                    metastatic.number), 
      "Primary Tumor"       = paste("Primary Tumor \n n= ", 
                                    primary.tumor.number),
      "Recurrent Tumor"     = paste("Recurrent Tumor \n n= ", 
                                    recurrent.tumor.number),
      "Solid Tissue Normal" = paste("Solid Tissue Normal \n n= ", 
                                    solid.tissue.normal.number))) +
    theme_bw() +
    theme(plot.title      = black_bold_tahoma_12,
          axis.title      = black_bold_tahoma_12,
          axis.text.x     = black_bold_tahoma_12_45,
          axis.text.y     = black_bold_tahoma_12,
          axis.line.x     = element_line(color = "black"),
          axis.line.y     = element_line(color = "black"),
          panel.grid      = element_blank(),
          legend.position = "none",
          strip.text      = black_bold_tahoma_12) +
    stat_compare_means(comparisons = list(
      c("Metastatic", "Solid Tissue Normal"), 
      c("Primary Tumor", "Solid Tissue Normal"), 
      c("Recurrent Tumor", "Solid Tissue Normal"),
      c("Metastatic", "Primary Tumor"),
      c("Recurrent Tumor", "Primary Tumor")), 
      method = "t.test")  
  
  p2 <- ggplot(data = x,
               aes(x     = variable, 
                   y     = value, 
                   color = variable)) +
    facet_grid( ~ sample.type, 
                scales = "free", 
                space  = "free")+   
    facet_wrap( ~ sample.type, ncol = 6)+
    geom_violin(trim = FALSE) +
    geom_boxplot(alpha    = .01, 
                 size     = .75,
                 width    = .5,
                 position = position_dodge(width = .9)) +
    labs(x = "EIF complex",
         y = paste("log2(value)")) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    theme_bw() +
    theme(plot.title      = black_bold_tahoma_12,
          axis.title      = black_bold_tahoma_12,
          axis.text.x     = black_bold_tahoma_12_45,
          axis.text.y     = black_bold_tahoma_12,
          axis.line.x     = element_line(color = "black"),
          axis.line.y     = element_line(color = "black"),
          panel.grid      = element_blank(),
          legend.position = "none",
          strip.text      = black_bold_tahoma_12) +
    stat_compare_means(comparisons = list(c("EIF4A1:EIF4E ratio", 
                                            "EIF4E:EIF4E ratio"), 
                                          c("EIF4G1:EIF4E ratio", 
                                            "EIF4E:EIF4E ratio"), 
                                          c("EIF4EBP1:EIF4E ratio", 
                                            "EIF4E:EIF4E ratio"),
                                          c("RPS6KB1:EIF4E ratio", 
                                            "EIF4E:EIF4E ratio"),
                                          c("MYC:EIF4E ratio", 
                                            "EIF4E:EIF4E ratio")), 
                       method = "t.test")  
print(p1)
print(p2)
  }

##
plotEIF.GTEX <-  function (x) {
  name <- deparse(substitute(x))
  cell.line.number <- nrow(x[x$sample.type == "Cell Line",])
  normal.tissue.number <- nrow(x[x$sample.type == "Normal Tissue",])
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
    facet_grid( ~ variable, 
                scales = "free", 
                space  = "free")+   
    facet_wrap( ~ variable, ncol = 6)+
    geom_violin(trim = FALSE) +
    geom_boxplot(alpha    = .01, 
                 size     = .75,
                 width    = .5,
                 position = position_dodge(width = .9)) +
    labs(x = "sample type",
         y = paste("log2(TPM)")) +
    scale_x_discrete(labels = c("Cell Line"     = paste("Cell Line \n n= ",
                                                        cell.line.number), 
                                "Normal Tissue" = paste("Normal Tissue \n n= ", 
                                                        normal.tissue.number))) +
    theme_bw() +
    theme(plot.title      = black_bold_tahoma_12,
          axis.title      = black_bold_tahoma_12,
          axis.text.x     = black_bold_tahoma_12_45,
          axis.text.y     = black_bold_tahoma_12,
          axis.line.x     = element_line(color = "black"),
          axis.line.y     = element_line(color = "black"),
          panel.grid      = element_blank(),
          legend.position = "none",
          strip.text      = black_bold_tahoma_12)+
    stat_compare_means(method = "anova")
}

##
plot.EIF.seq.all.samples <- function(x){
  plotEIF.TCGA.GTEX(x)+
    stat_compare_means(method = "anova")
} 

##
plot.EIF.seq.each.tumor <- function(x, y){
  m <- x[x$primary.disease == y,]
  my_comparison <- list(c("Metastatic", "Solid Tissue Normal"), 
                        c("Primary Tumor", "Solid Tissue Normal"), 
                        c("Recurrent Tumor", "Solid Tissue Normal"),
                        c("Metastatic", "Primary Tumor"),
                        c("Recurrent Tumor", "Primary Tumor"))
  plotEIF.RNAseq.TCGA(m)+
    labs(title = y) +
    stat_compare_means(method = "anova")
} 

##
get.disease.list <- function () {
  x <- get.EIF.TCGA.RNAseq.long()
  disease.list <- levels(x$primary.disease)
  names(disease.list)<- disease.list
  return(disease.list)}

##
plot.EIF.seq.all.normal <- function(x){
  plotEIF.GTEX(x)+
    stat_compare_means(method = "t.test")
} 

##
plot.EIF.PCA <- function (){
  EIF.proteomics <- read.csv(file.path("project-data", 
                                       "proteomics.csv"), 
                             header = TRUE, 
                             sep = ",")
  EIF.proteomics <- as.data.frame(EIF.proteomics)
#  EIF.proteomics$gene_id <- as.factor(EIF.proteomics$gene_id)
  EIF.proteomics <- na.omit(EIF.proteomics)
  df1 <- EIF.proteomics[,-1]
  df <- as.numeric(df)
  ggplot2::autoplot(prcomp(df1))

  ggplot2::autoplot(prcomp(df), 
                    data                 = DNFA.TCGA.GTEX.skin, 
                    colour               = 'sample_type',
                    loadings             = TRUE, 
                    loadings.colour      = 'black',
                    loadings.label.vjust = 0,
                    loadings.label       = TRUE, 
                    loadings.label.size  = 6) +
    theme(plot.background  = element_blank(),
          panel.background = 
            element_rect(fill  = 'transparent',
                         color = 'black',
                         size  = 1),
          axis.title   = element_text(colour  = "black", 
                                      size    = 24, 
                                      face    = "bold"),
          axis.text    = element_text(colour  = "black", 
                                      size    = 24, 
                                      face    = "bold"),
          legend.title = element_text(colour  = "black", 
                                      size    = 24, 
                                      face    = "bold"),
          legend.text  = element_text(colour  = "black", 
                                      size    = 24, 
                                      face    = "bold",
                                      hjust   = 1),
          legend.key   = element_blank())
}

##
plot.DNFA.PCA <- function (){
  DNFA.TCGA.GTEX <- read.csv(file.path("project-data", 
                                       "DNFATCGAandGTEX.csv"), 
                             header = TRUE, 
                             sep = ",")
  DNFA.TCGA.GTEX <- as.data.frame(DNFA.TCGA.GTEX)
  DNFA.TCGA.GTEX$sample_type <- as.factor(DNFA.TCGA.GTEX$sample_type)
  DNFA.TCGA.GTEX$sample_type <- factor(DNFA.TCGA.GTEX$sample_type, 
                                         levels = c("Normal Tissue", 
                                                    "Primary Tumor", 
                                                    "Metastatic", 
                                                    "Solid Tissue Normal"))
  sample.type <- c("Skin")
  DNFA.TCGA.GTEX.skin <- DNFA.TCGA.GTEX[
    DNFA.TCGA.GTEX$primary.disease.or.tissue %in% sample.type, ]
  DNFA.TCGA.GTEX.skin <- na.omit(DNFA.TCGA.GTEX.skin)
  df <- DNFA.TCGA.GTEX.skin[c(3,4,5,6,7,9)]
  df <- na.omit(df)
  ggplot2::autoplot(prcomp(df))
  # autoplot(prcomp(df), data = DNFASKCMandTGEX, colour = 'sample_type')
  
  ggplot2::autoplot(prcomp(df), 
                    data                 = DNFA.TCGA.GTEX.skin, 
                    colour               = 'sample_type',
                    loadings             = TRUE, 
                    loadings.colour      = 'black',
                    loadings.label.vjust = 0,
                    loadings.label       = TRUE, 
                    loadings.label.size  = 6) +
    theme(plot.background  = element_blank(),
          panel.background = 
            element_rect(fill  = 'transparent',
                         color = 'black',
                         size  = 1),
          axis.title   = element_text(colour  = "black", 
                                      size    = 24, 
                                      face    = "bold"),
          axis.text    = element_text(colour  = "black", 
                                      size    = 24, 
                                      face    = "bold"),
          legend.title = element_text(colour  = "black", 
                                      size    = 24, 
                                      face    = "bold"),
          legend.text  = element_text(colour  = "black", 
                                      size    = 24, 
                                      face    = "bold",
                                      hjust   = 1),
          legend.key   = element_blank())
}
####################################################################
##  Kaplan-Meier curve with clinic and EIF RNASeq data all tumor  ##
####################################################################
plot.km.EIF.all.tumors <- function(EIF) {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data", "EIFTCGAGTEX2.csv"), 
                            header = TRUE, sep = ",")
  EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$study == 'TCGA',]
  EIF.TCGA <- EIF.TCGA[EIF.TCGA$X_sample_type != "Solid Tissue Normal", ]
  EIF.TCGA <- droplevels(EIF.TCGA)
  df <- na.omit(EIF.TCGA)
  number <- nrow(df)
  sub <- round(number/5, digits = 0)
  bottom.label <- paste("Bottom 20%, n = ", sub)
  top.label <- paste("Top 20%, n = ", sub)
  df$Group[df[[EIF]] < quantile(df[[EIF]], prob = 0.2)] = "Bottom 20%"
  df$Group[df[[EIF]] > quantile(df[[EIF]], prob = 0.8)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS.time, OS == 1))
  df <- na.omit(df)
  km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 0) # rho = 0 log-rank
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
  black.bold.12pt <- element_text(face   = "bold",
                                  size   = 12,
                                  family = "Tahoma", 
                                  colour = "black")
  print(
    ggplot2::autoplot(km,
                      xlim = c(0, 4000),
                      xlab = "Days",
                      ylab = "Survival Probability",
                      main = paste0("Kaplan-Meier plot of all TCGA cancer studies(", 
                                    number," cases)")) +
      theme_bw() +
      theme(plot.title           = black.bold.12pt,
            axis.title           = black.bold.12pt,
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
                         labels = c(bottom.label, top.label)) +
      geom_point(size = 0.25) +
      annotate("text",
               x        = 4000,
               y        = 0.8,
               label    = paste("log-rank test, p.val = ", p.val),
               size     = 4.5,
               hjust    = 1,
               fontface = "bold"))
  # rho = 1 the Gehan-Wilcoxon test
  print(EIF)
  print(stats)
  #  fit = survfit(SurvObj ~ df$Group, data = df)
  #  tst <- comp(fit)$tests$lrTests
  #  print(tst)
}

##
plot.km.DNFA.all.tumors <- function(DNFA) {
  DNFA.TCGA.GTEX <- read.csv(file.path("project-data", "DNFASKCMandGTEX.csv"), 
                            header = TRUE, sep = ",")
  DNFA.TCGA <- DNFA.TCGA.GTEX[DNFA.TCGA.GTEX$study == 'TCGA',]
  DNFA.TCGA <- DNFA.TCGA[DNFA.TCGA$sample_type != "Solid Tissue Normal", ]
  DNFA.TCGA <- droplevels(DNFA.TCGA)
  df <- na.omit(DNFA.TCGA)
  number <- nrow(df)
  sub <- round(number/5, digits = 0)
  bottom.label <- paste("Bottom 10%, n = ", sub)
  top.label <- paste("Top 10%, n = ", sub)
  df$Group[df[[DNFA]] < quantile(df[[DNFA]], prob = 0.1)] = "Bottom"
  df$Group[df[[DNFA]] > quantile(df[[DNFA]], prob = 0.9)] = "Top"
  df$SurvObj <- with(df, Surv(OS.time, OS == 1))
  df <- na.omit(df)
  km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 0) # rho = 0 log-rank
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
  black.bold.12pt <- element_text(face   = "bold",
                                  size   = 12,
                                  family = "Tahoma", 
                                  colour = "black")
  print(
    ggplot2::autoplot(km,
                      xlab = "Days",
                      xlim = c(0, 4000),
                      ylab = "Survival Probability",
                      main = paste0("Kaplan-Meier plot of all TCGA cancer studies(", 
                                    number," cases)")) +
      theme_bw() +
      theme(plot.title           = black.bold.12pt,
            axis.title           = black.bold.12pt,
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
                         name   = paste(DNFA, "mRNA expression"),
                         breaks = c("Bottom", "Top"),
                         labels = c(bottom.label, top.label)) +
      geom_point(size = 0.25) +
      annotate("text",
               x        = 4000,
               y        = 0.8,
               label    = paste("log-rank test, p.val = ", p.val),
               size     = 4.5,
               hjust    = 1,
               fontface = "bold"))
  # rho = 1 the Gehan-Wilcoxon test
  print(DNFA)
  print(stats)
  #  fit = survfit(SurvObj ~ df$Group, data = df)
  #  tst <- comp(fit)$tests$lrTests
  #  print(tst)
}

###################################################################################
##  Kaplan-Meier curve with clinic and EIF RNASeq data in individual tumor group ##
###################################################################################
plot.km.EIF.each.tumor <- function(EIF, tumor) {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data", "EIFTCGAGTEX2.csv"), 
                            header = TRUE, sep = ",")
  EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$study == 'TCGA',]
  EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$primary.disease.or.tissue == tumor,]
  EIF.TCGA <- EIF.TCGA[EIF.TCGA$X_sample_type != "Solid Tissue Normal", ]
  EIF.TCGA <- droplevels(EIF.TCGA)
  df <- na.omit(EIF.TCGA)
  number <- nrow(df)
  sub <- round(number/5, digits = 0)
  bottom.label <- paste("Bottom 20%, n = ", sub)
  top.label <- paste("Top 20%, n = ", sub)
  df$Group[df[[EIF]] < quantile(df[[EIF]], prob = 0.2)] = "Bottom 20%"
  df$Group[df[[EIF]] > quantile(df[[EIF]], prob = 0.8)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS.time, OS == 1))
  df <- na.omit(df)
  km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 0) # rho = 0 log-rank
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
  black.bold.12pt <- element_text(face   = "bold",
                                  size   = 12,
                                  family = "Tahoma", 
                                  colour = "black")
  print(
    ggplot2::autoplot(km,
                      xlim = c(0, 4000),
                      xlab = "Days",
                      ylab = "Survival Probability",
                      main = paste0("Kaplan-Meier plot of ", tumor, " (",
                                    number," cases)")) +
      theme_bw() +
      theme(plot.title           = black.bold.12pt,
            axis.title           = black.bold.12pt,
            axis.text            = black.bold.12pt,
            axis.line.x          = element_line(color = "black"),
            axis.line.y          = element_line(color = "black"),
            panel.grid           = element_blank(),
            strip.text           = black.bold.12pt,
            legend.text          = black.bold.12pt ,
            legend.title         = black.bold.12pt ,
            legend.position      = c(1,1),
            legend.justification = c(1,1)) +
      guides(fill = FALSE) +
      scale_color_manual(values  = c("red", "blue"),
                         name    = paste(EIF, "mRNA expression"),
                         breaks  = c("Bottom 20%", "Top 20%"),
                         labels  = c(bottom.label, top.label)) +
      geom_point(size = 0.25) +
      annotate("text",
               x        = 4000,
               y        = 0.8,
               label    = paste("log-rank test, p.val = ", p.val),
               size     = 4.5,
               hjust    = 1,
               fontface = "bold"))
  # rho = 1 the Gehan-Wilcoxon test
  print(EIF)
  print(stats)
  #  fit = survfit(SurvObj ~ df$Group, data = df)
  #  tst <- comp(fit)$tests$lrTests
  #  print(tst)
}

##
plot.km.DNFA.each.tumor <- function(EIF, tumor) {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data", "DNFASKCMandGTEX.csv"), 
                            header = TRUE, sep = ",")
  EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$study == 'TCGA',]
  EIF.TCGA <- EIF.TCGA[EIF.TCGA$primary.disease.or.tissue == tumor,]
  EIF.TCGA <- EIF.TCGA[EIF.TCGA$sample_type != "Solid Tissue Normal", ]
#  EIF.TCGA <- EIF.TCGA[which(EIF.TCGA$OS.time < 4001), ]
#  EIF.TCGA <- subset(EIF.TCGA,  OS.time < 4000)
#  EIF.TCGA <- subset(EIF.TCGA,  OS.time > 0)
  EIF.TCGA <- droplevels(EIF.TCGA)
  df <- na.omit(EIF.TCGA)
  number <- nrow(df)
  sub <- round(number/5, digits = 0)
  bottom.label <- paste("Bottom 30%, n = ", sub)
  top.label <- paste("Top 30%, n = ", sub)
  df$Group[df[[EIF]] < quantile(df[[EIF]], prob = 0.3)] = "Bottom 20%"
  df$Group[df[[EIF]] > quantile(df[[EIF]], prob = 0.7)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS.time, OS == 1))
  df <- na.omit(df)
  km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 0) # rho = 0 log-rank
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
#  p.val.2 <- surv_pvalue(km)
  black.bold.12pt <- element_text(face   = "bold",
                                  size   = 12,
                                  family = "Tahoma", 
                                  colour = "black")
  print(
    ggplot2::autoplot(km,
                      xlab = "Days",
                      xlim = c(0, 4000),
                      ylab = "Survival Probability",
                      main = paste0("Kaplan-Meier plot of ", tumor, " (",
                                    number," cases)")) +
      theme_bw() +
      theme(plot.title           = black.bold.12pt,
            axis.title           = black.bold.12pt,
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
                         labels = c(bottom.label, top.label)) +
      geom_point(size = 0.25) +
      annotate("text",
               x        = 4000,
               y        = 0.8,
               label    = paste("log-rank test, p.val = ", p.val),
               size     = 4.5,
               hjust    = 1,
               fontface = "bold"))
  # rho = 1 the Gehan-Wilcoxon test
  print(EIF)
  print(stats)
  #  fit = survfit(SurvObj ~ df$Group, data = df)
  #  tst <- comp(fit)$tests$lrTests
  #  print(tst)
}

####################################################
####################################################
plotEIF.TCGA.GTEX (get.EIF.TCGA.GTEX.RNAseq.long())
plotEIF.TCGA.GTEX (get.EIF.TCGA.GTEX.score.long())

plotEIF.TCGA.GTEX (get.DNFA.TCGA.GTEX.RNAseq.long())

plotEIF.RNAseq.TCGA (get.EIF.TCGA.RNAseq.long())
plotEIF.score.TCGA (get.EIF.TCGA.score.long())

plot.EIF.seq.all.normal (get.EIF.GTEX.RNAseq.long())
plot.EIF.seq.all.normal (get.EIF.GTEX.score.long())

plot.EIF.seq.each.tumor (x = get.EIF.TCGA.RNAseq.long(), 
                         y = "Prostate Adenocarcinoma")
lapply(get.disease.list(), 
       plot.EIF.seq.each.tumor, 
       x = get.EIF.TCGA.RNAseq.long())

plot.EIF.seq.each.tumor (x = get.DNFA.TCGA.GTEX.RNAseq.long(), 
                         y = "Skin Cutaneous Melanoma")
lapply(get.disease.list(), 
       plot.EIF.seq.each.tumor, 
       x = get.DNFA.TCGA.GTEX.RNAseq.long())
####################################################
####################################################
plot.km.DNFA.all.tumors("HMGCR")
plot.km.EIF.all.tumors("SF3B5")

DNFA.gene <- c("ACLY", "ACSS2","ACACA", "SCD", "FASN", "ACSL1", 
               "HMGCS1", "HMGCR", "MVK", "PMVK")
names(DNFA.gene ) <- DNFA.gene

EIF.gene <- c("EIF4A1","EIF4E","EIF4G1","EIF4EBP1","RPS6KB1","MYC",
              "HIST2H2BE", "HIST1H1B", "HIST1H2AB","HIST1H1C",
              "PPA2","ECHS1", "LLGL2","SF3B5", "RPS5")
names(EIF.gene) <- EIF.gene



sapply(DNFA.gene, plot.km.DNFA.all.tumors)
sapply(EIF.gene, plot.km.EIF.all.tumors)

plot.km.EIF.each.tumor("EIF4EBP1", "Cervical & Endocervical Cancer")
lapply(get.disease.list(), 
       plot.km.EIF.each.tumor, 
       EIF = "EIF4A1")

lapply(EIF.gene, 
       plot.km.EIF.each.tumor, 
       tumor = "Cervical & Endocervical Cancer")


lapply(get.disease.list(), 
       plot.km.DNFA.each.tumor, 
       EIF = "SCD")


plot.km.DNFA.each.tumor("HMGCR", "Skin Cutaneous Melanoma")
lapply(DNFA.gene, 
       plot.km.DNFA.each.tumor, 
       tumor = "Bladder Urothelial Carcinoma")

lapply(DNFA.gene, 
       plot.km.DNFA.each.tumor, 
       tumor = "Sarcoma")

lapply(DNFA.gene, 
       plot.km.DNFA.each.tumor, 
       tumor = "Cervical & Endocervical Cancer")




















