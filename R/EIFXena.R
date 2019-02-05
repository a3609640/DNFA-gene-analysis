## the following script perform PCA on RNA-Seq data of seven DNFA genes 
## from SKCM amd GTEX dataset with R package "ggfortify".

library(ggfortify)
library(ggplot2)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(gridExtra)
library(reshape2)

## read.csv will transform characters into factors  
get.EIF.TCGA.GTEX <- function () {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data", "EIFTCGAGTEX.csv"), 
                            header = TRUE, sep = ",")
  EIF.TCGA.GTEX.RNAseq.long <- melt(EIF.TCGA.GTEX[, 1:10])
  colnames(EIF.TCGA.GTEX.RNAseq.long) <- c("sample", "study", "sample.type", 
                                           "primary.disease", "variable", "value")
  return(EIF.TCGA.GTEX.RNAseq.long)
}

##
get.EIF.TCGA.RNAseq.long <- function () {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data", 
                                      "EIFTCGAGTEX.csv"), 
                            header = TRUE, sep = ",")
  EIF.TCGA.GTEX.RNAseq.long <- melt(EIF.TCGA.GTEX[, 1:10])
  colnames(EIF.TCGA.GTEX.RNAseq.long) <- c("sample", "study", "sample.type", 
                                           "primary.disease", "variable", "value")
  EIF.TCGA.RNAseq.long <- EIF.TCGA.GTEX.RNAseq.long[EIF.TCGA.GTEX.RNAseq.long$study == 'TCGA',]
  EIF.TCGA.RNAseq.long <- EIF.TCGA.RNAseq.long[EIF.TCGA.RNAseq.long$value != 0,]
  tumor.type <- c("Metastatic","Primary Tumor",
                  "Recurrent Tumor","Normal Tissue",
                  "Solid Tissue Normal")
  EIF.TCGA.RNAseq.long <- EIF.TCGA.RNAseq.long[EIF.TCGA.RNAseq.long$sample.type %in% tumor.type,]
  EIF.TCGA.RNAseq.long <- droplevels(EIF.TCGA.RNAseq.long)
  return(EIF.TCGA.RNAseq.long)
}

##
plotEIF <-  function (x) {
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
  ggplot(data = x,
         aes(x     = sample.type, 
             y     = value, 
             color = sample.type)) +
    facet_grid( ~ variable, 
                scales = "free", 
                space  = "free")+   
    facet_wrap( ~ variable, ncol = 3)+
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
          strip.text      = black_bold_tahoma_12)
}

##
plot.EIF.seq.all.tumors <- function(x){
  my_comparison <- list(c("Metastatic", "Solid Tissue Normal"), 
                        c("Primary Tumor", "Solid Tissue Normal"), 
                        c("Recurrent Tumor", "Solid Tissue Normal"),
                        c("Metastatic", "Primary Tumor"),
                        c("Recurrent Tumor", "Primary Tumor"))
  plotEIF(x)+
    stat_compare_means(comparisons = my_comparison, method = "t.test")
} 

##
plot.EIF.seq.each.tumor <- function(x, y){
  m <- x[x$primary.disease == y,]
  my_comparison <- list(c("Metastatic", "Solid Tissue Normal"), 
                        c("Primary Tumor", "Solid Tissue Normal"), 
                        c("Recurrent Tumor", "Solid Tissue Normal"),
                        c("Metastatic", "Primary Tumor"),
                        c("Recurrent Tumor", "Primary Tumor"))
  plotEIF(m)+
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




####################################################################
##  Kaplan-Meier curve with clinic and EIF RNASeq data from hnsc  ##
####################################################################
plot.km.EIF.all.tumors <- function(EIF) {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data", "EIFTCGAGTEX.csv"), 
                            header = TRUE, sep = ",")
  EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$X_study == 'TCGA',]
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
                                  colour = "black")
  print(
    ggplot2::autoplot(km,
                      xlab = "Months",
                      ylab = "Survival Probability",
                      main = paste0("Kaplan-Meier plot (", 
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
                         labels = c(bottom.label, top.label)) +
      geom_point(size = 0.25) +
      annotate("text",
               x        = 10000,
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

plot.km.EIF.each.tumor <- function(EIF, tumor) {
  EIF.TCGA.GTEX <- read.csv(file.path("project-data", "EIFTCGAGTEX.csv"), 
                            header = TRUE, sep = ",")
  EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$X_study == 'TCGA',]
  EIF.TCGA <- EIF.TCGA.GTEX[EIF.TCGA.GTEX$primary.disease.or.tissue == tumor,]
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
                                  colour = "black")
  print(
    ggplot2::autoplot(km,
                      xlab = "Days",
                      ylab = "Survival Probability",
                      main = paste0("Kaplan-Meier plot of ", tumor, " (",
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
                         labels = c(bottom.label, top.label)) +
      geom_point(size = 0.25) +
      annotate("text",
               x        = 10000,
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

plot.km.EIF.all.tumors("EIF4G1")
EIF.gene <- c("EIF4A1","EIF4E","EIF4G1","EIF4EBP1","RPS6KB1","MYC")
names(EIF.gene) <- EIF.gene
lapply(EIF.gene, plot.km.EIF.all.tumors)



plot.km.EIF.each.tumor ("EIF4G1", "Sarcoma")
lapply(get.disease.list(), plot.km.EIF.each.tumor, EIF = "EIF4G1")















####################################################



plot.EIF.seq.all.tumors (get.EIF.TCGA.RNAseq.long())
plot.EIF.seq.each.tumor (x = get.EIF.TCGA.RNAseq.long(), y = "Sarcoma")
lapply(get.disease.list(), plot.each.tumor.type, x = get.EIF.TCGA.RNAseq.long())













