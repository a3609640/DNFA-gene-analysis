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



####################################################
## Boxplots of relative levels of DNFA expression ##
####################################################
## read.csv will transform characters into factors
get.DNFA.TCGA.GTEX.RNAseq.long <- function () {
  DNFA.TCGA.GTEX <- read.csv(
    file.path("project-data",
      "DNFASKCMandGTEX.csv"),
    header = TRUE,
    sep = ","
  )
  DNFA.TCGA.GTEX.RNAseq.long <- melt(DNFA.TCGA.GTEX[, 1:14])
  colnames(DNFA.TCGA.GTEX.RNAseq.long) <- c("sample",
    "study",
    "sample.type",
    "primary.disease",
    "variable",
    "value")
  DNFA.TCGA.GTEX.RNAseq.long <- DNFA.TCGA.GTEX.RNAseq.long[DNFA.TCGA.GTEX.RNAseq.long$value != 0, ]
  DNFA.TCGA.GTEX.RNAseq.long <- na.omit(DNFA.TCGA.GTEX.RNAseq.long)
  tumor.type <- c(
    "Metastatic",
    "Primary Tumor",
    "Recurrent Tumor",
    "Solid Tissue Normal",
    "Normal Tissue",
    "Cell Line"
  )
  DNFA.TCGA.GTEX.RNAseq.long <- DNFA.TCGA.GTEX.RNAseq.long[DNFA.TCGA.GTEX.RNAseq.long$sample.type %in% tumor.type, ]
  DNFA.TCGA.GTEX.RNAseq.long <-
    droplevels(DNFA.TCGA.GTEX.RNAseq.long)
  DNFA.TCGA.GTEX.RNAseq.long$sample.type <-
    factor(DNFA.TCGA.GTEX.RNAseq.long$sample.type,
      levels = tumor.type)
  return(DNFA.TCGA.GTEX.RNAseq.long)
}

get.DNFA.TCGA.RNAseq.long <- function () {
  DNFA.TCGA.GTEX <- read.csv(
    file.path("project-data",
      "DNFASKCMandGTEX.csv"),
    header = TRUE,
    sep = ","
  )
  DNFA.TCGA.GTEX.RNAseq.long <- melt(DNFA.TCGA.GTEX[, 1:14])
  colnames(DNFA.TCGA.GTEX.RNAseq.long) <- c("sample",
    "study",
    "sample.type",
    "primary.disease",
    "variable",
    "value")
  DNFA.TCGA.RNAseq.long <- DNFA.TCGA.GTEX.RNAseq.long[DNFA.TCGA.GTEX.RNAseq.long$study == 'TCGA', ]
  DNFA.TCGA.RNAseq.long <- DNFA.TCGA.RNAseq.long[DNFA.TCGA.RNAseq.long$value != 0, ]
  tumor.type <- c(
    "Metastatic",
    "Primary Tumor",
    "Recurrent Tumor",
    "Normal Tissue",
    "Solid Tissue Normal"
  )
  DNFA.TCGA.RNAseq.long <- DNFA.TCGA.RNAseq.long[DNFA.TCGA.RNAseq.long$sample.type %in% tumor.type, ]
  DNFA.TCGA.RNAseq.long <- droplevels(DNFA.TCGA.RNAseq.long)
  return(DNFA.TCGA.RNAseq.long)
}

get.DNFA.GTEX.RNAseq.long <- function () {
  DNFA.TCGA.GTEX <- read.csv(
    file.path("project-data",
      "DNFASKCMandGTEX.csv"),
    header = TRUE,
    sep = ","
  )
  DNFA.TCGA.GTEX.RNAseq.long <- melt(DNFA.TCGA.GTEX[, 1:14])
  colnames(DNFA.TCGA.GTEX.RNAseq.long) <- c("sample",
    "study",
    "sample.type",
    "primary.disease",
    "variable",
    "value")
  DNFA.GTEX.RNAseq.long <- DNFA.TCGA.GTEX.RNAseq.long[DNFA.TCGA.GTEX.RNAseq.long$study == 'GTEX', ]
  DNFA.GTEX.RNAseq.long <- DNFA.GTEX.RNAseq.long[DNFA.GTEX.RNAseq.long$value != 0, ]
  DNFA.GTEX.RNAseq.long <- na.omit(DNFA.GTEX.RNAseq.long)
  #  DNFA.GTEX.RNAseq.long <- DNFA.GTEX.RNAseq.long[DNFA.GTEX.RNAseq.long$sample.type %in% tumor.type,]
  DNFA.GTEX.RNAseq.long <- droplevels(DNFA.GTEX.RNAseq.long)
  return(DNFA.GTEX.RNAseq.long)
}

### plot function for get.DNFA.TCGA.GTEX.RNAseq.long()
plot.DNFA.TCGA.GTEX <-  function (x) {
  name <- deparse(substitute(x))
  metastatic.number <- nrow(x[x$sample.type == "Metastatic", ])
  primary.tumor.number <- nrow(x[x$sample.type == "Primary Tumor", ])
  recurrent.tumor.number <-
    nrow(x[x$sample.type == "Recurrent Tumor", ])
  solid.tissue.normal.number <-
    nrow(x[x$sample.type == "Solid Tissue Normal", ])
  cell.line.number <- nrow(x[x$sample.type == "Cell Line", ])
  normal.tissue.number <- nrow(x[x$sample.type == "Normal Tissue", ])
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 9
  )
  black_bold_tahoma_12_45 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 9,
    angle  = 45,
    hjust  = 1
  )
  p1 <- ggplot(data = x,
    aes(x     = sample.type,
      y     = value,
      color = sample.type)) +
    facet_grid(~ variable,
      scales = "free",
      space  = "free") +
    facet_wrap(~ variable, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = "sample type",
      y = paste("log2(TPM)")) +
    scale_x_discrete(
      labels =
        c(
          "Metastatic"          = paste("Metastatic \n n= ",
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
            cell.line.number)
        )
    ) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    )
  p1 <- p1 + stat_compare_means(method = "anova")
  p2 <- ggplot(data = x,
    aes(x     = variable,
      y     = value,
      color = variable)) +
    facet_grid(~ sample.type,
      scales = "free",
      space  = "free") +
    facet_wrap(~ sample.type, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = "sample type",
      y = paste("log2(value)")) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    )
  p2 <- p2 + stat_compare_means(method = "anova")
  print(p1)
  print(p2)
}

### plot function for get.DNFA.TCGA.RNAseq.long()
plot.DNFA.TCGA <-  function (x) {
  name <- deparse(substitute(x))
  metastatic.number <- nrow(x[x$sample.type == "Metastatic", ])
  primary.tumor.number <- nrow(x[x$sample.type == "Primary Tumor", ])
  recurrent.tumor.number <-
    nrow(x[x$sample.type == "Recurrent Tumor", ])
  solid.tissue.normal.number <-
    nrow(x[x$sample.type == "Solid Tissue Normal", ])
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12
  )
  black_bold_tahoma_12_45 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12,
    angle  = 45,
    hjust  = 1
  )
  p1 <- ggplot(data = x,
    aes(x     = sample.type,
      y     = value,
      color = sample.type)) +
    facet_grid(~ variable,
      scales = "free",
      space  = "free") +
    facet_wrap(~ variable, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = "sample type",
      y = paste("log2(RNA counts)")) +
    scale_x_discrete(
      labels =
        c(
          "Metastatic"          = paste("Metastatic \n n= ",
            metastatic.number),
          "Primary Tumor"       = paste("Primary Tumor \n n= ",
            primary.tumor.number),
          "Recurrent Tumor"     = paste("Recurrent Tumor \n n= ",
            recurrent.tumor.number),
          "Solid Tissue Normal" = paste("Solid Tissue Normal \n n= ",
            solid.tissue.normal.number)
        )
    ) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    ) +
    stat_compare_means(comparisons = list(
      c("Metastatic",
        "Solid Tissue Normal"),
      c("Primary Tumor",
        "Solid Tissue Normal"),
      c("Recurrent Tumor",
        "Solid Tissue Normal"),
      c("Metastatic",
        "Primary Tumor"),
      c("Recurrent Tumor",
        "Primary Tumor")
    ),
      method = "t.test")
  
  p2 <- ggplot(data = x,
    aes(x     = variable,
      y     = value,
      color = variable)) +
    facet_grid(~ sample.type,
      scales = "free",
      space  = "free") +
    facet_wrap(~ sample.type, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = "EIF complex",
      y = paste("log2(RNA counts)")) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    ) +
    stat_compare_means(method = "anova")
  print(p1)
  print(p2)
}

### plot function for get.DNFA.GTEX.RNAseq.long()
plot.DNFA.GTEX <-  function (x) {
  name <- deparse(substitute(x))
  cell.line.number <- nrow(x[x$sample.type == "Cell Line", ])
  normal.tissue.number <- nrow(x[x$sample.type == "Normal Tissue", ])
  black_bold_tahoma_12 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12
  )
  black_bold_tahoma_12_45 <- element_text(
    color  = "black",
    face   = "bold",
    family = "Tahoma",
    size   = 12,
    angle  = 45,
    hjust  = 1
  )
  ggplot(data = x,
    aes(x     = sample.type,
      y     = value,
      color = sample.type)) +
    facet_grid(~ variable,
      scales = "free",
      space  = "free") +
    facet_wrap(~ variable, ncol = 6) +
    geom_violin(trim = FALSE) +
    geom_boxplot(
      alpha    = .01,
      size     = .75,
      width    = .5,
      position = position_dodge(width = .9)
    ) +
    labs(x = "sample type",
      y = paste("log2(TPM)")) +
    scale_x_discrete(labels = c(
      "Cell Line"     = paste("Cell Line \n n= ",
        cell.line.number),
      "Normal Tissue" = paste("Normal Tissue \n n= ",
        normal.tissue.number)
    )) +
    theme_bw() +
    theme(
      plot.title      = black_bold_tahoma_12,
      axis.title      = black_bold_tahoma_12,
      axis.text.x     = black_bold_tahoma_12_45,
      axis.text.y     = black_bold_tahoma_12,
      axis.line.x     = element_line(color = "black"),
      axis.line.y     = element_line(color = "black"),
      panel.grid      = element_blank(),
      legend.position = "none",
      strip.text      = black_bold_tahoma_12
    ) +
    stat_compare_means(method = "anova")
}

### plot function for get.DNFA.TCGA.RNAseq.long(),
### allow selection of individual tumor group
plot.DNFA.TCGA.each.tumor <- function(x, y) {
  m <- x[x$primary.disease == y, ]
  my_comparison <- list(
    c("Metastatic", "Solid Tissue Normal"),
    c("Primary Tumor", "Solid Tissue Normal"),
    c("Recurrent Tumor", "Solid Tissue Normal"),
    c("Metastatic", "Primary Tumor"),
    c("Recurrent Tumor", "Primary Tumor")
  )
  plot.DNFA.TCGA(m) +
    labs(title = y) +
    stat_compare_means(method = "anova")
}

get.disease.list <- function () {
  x <- get.DNFA.TCGA.RNAseq.long()
  disease.list <- levels(x$primary.disease)
  names(disease.list) <- disease.list
  return(disease.list)
}

###########################################################################
## PCA with DNFA or DNCS gene expression on SKCM and normal skin tissues ##
###########################################################################
plot.DNFA.PCA <- function (x) {
  DNFA.TCGA.GTEX <- read.csv(
    file.path("project-data",
      "DNFATCGAandGTEX.csv"),
    header = TRUE,
    sep = ","
  )
  DNFA.TCGA.GTEX <- as.data.frame(DNFA.TCGA.GTEX)
  DNFA.TCGA.GTEX$sample_type <-
    as.factor(DNFA.TCGA.GTEX$sample_type)
  DNFA.TCGA.GTEX$sample_type <- factor(
    DNFA.TCGA.GTEX$sample_type,
    levels = c(
      "Normal Tissue",
      "Primary Tumor",
      "Metastatic",
      "Solid Tissue Normal"
    )
  )
  sample.type <- c("Skin")
  DNFA.TCGA.GTEX.skin <- DNFA.TCGA.GTEX[DNFA.TCGA.GTEX$primary.disease.or.tissue %in% sample.type,]
  DNFA.TCGA.GTEX.skin <- na.omit(DNFA.TCGA.GTEX.skin)
  df <- DNFA.TCGA.GTEX.skin[x]
  df <- na.omit(df)
  ggplot2::autoplot(prcomp(df))
  # autoplot(prcomp(df), data = DNFASKCMandTGEX, colour = 'sample_type')
  
  p <- ggplot2::autoplot(
    prcomp(df),
    data                 = DNFA.TCGA.GTEX.skin,
    colour               = 'sample_type',
    loadings             = TRUE,
    loadings.colour      = 'black',
    loadings.label.vjust = 0,
    loadings.label       = TRUE,
    loadings.label.size  = 6
  ) +
    theme(
      plot.background  = element_blank(),
      panel.background =
        element_rect(
          fill  = 'transparent',
          color = 'black',
          size  = 1
        ),
      axis.title   = element_text(
        colour  = "black",
        size    = 24,
        face    = "bold"
      ),
      axis.text    = element_text(
        colour  = "black",
        size    = 24,
        face    = "bold"
      ),
      legend.title = element_text(
        colour  = "black",
        size    = 24,
        face    = "bold"
      ),
      legend.text  = element_text(
        colour  = "black",
        size    = 24,
        face    = "bold",
        hjust   = 1
      ),
      legend.key   = element_blank()
    )
  print(p)
}

####################################################################
##  Kaplan-Meier curve with clinic and DNFA RNASeq data all tumor ##
####################################################################
plot.km.DNFA.all.tumors <- function(DNFA) {
  DNFA.TCGA.GTEX <-
    read.csv(
      file.path("project-data", "DNFASKCMandGTEX.csv"),
      header = TRUE,
      sep = ","
    )
  DNFA.TCGA <- DNFA.TCGA.GTEX[DNFA.TCGA.GTEX$study == 'TCGA', ]
  DNFA.TCGA <-
    DNFA.TCGA[DNFA.TCGA$sample_type != "Solid Tissue Normal",]
  DNFA.TCGA <- droplevels(DNFA.TCGA)
  df <- na.omit(DNFA.TCGA)
  number <- nrow(df)
  sub <- round(number / 5, digits = 0)
  bottom.label <- paste("Bottom 20%, n = ", sub)
  top.label <- paste("Top 20%, n = ", sub)
  df$Group[df[[DNFA]] < quantile(df[[DNFA]], prob = 0.2)] = "Bottom"
  df$Group[df[[DNFA]] > quantile(df[[DNFA]], prob = 0.8)] = "Top"
  df$SurvObj <- with(df, Surv(OS.time, OS == 1))
  df <- na.omit(df)
  km <-
    survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <-
    survdiff(SurvObj ~ df$Group, data = df, rho = 0) # rho = 0 log-rank
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
  black.bold.12pt <- element_text(
    face   = "bold",
    size   = 12,
    family = "Tahoma",
    colour = "black"
  )
  print(
    ggplot2::autoplot(
      km,
      xlab = "Days",
      xlim = c(0, 4000),
      ylab = "Survival Probability",
      main = paste0(
        "Kaplan-Meier plot of all TCGA cancer studies(",
        number,
        " cases)"
      )
    ) +
      theme_bw() +
      theme(
        plot.title           = black.bold.12pt,
        axis.title           = black.bold.12pt,
        axis.text            = black.bold.12pt,
        axis.line.x          = element_line(color  = "black"),
        axis.line.y          = element_line(color  = "black"),
        panel.grid           = element_blank(),
        strip.text           = black.bold.12pt,
        legend.text          = black.bold.12pt ,
        legend.title         = black.bold.12pt ,
        legend.position      = c(1, 1),
        legend.justification = c(1, 1)
      ) +
      guides(fill = FALSE) +
      scale_color_manual(
        values = c("red", "blue"),
        name   = paste(DNFA, "mRNA expression"),
        breaks = c("Bottom", "Top"),
        labels = c(bottom.label, top.label)
      ) +
      geom_point(size = 0.25) +
      annotate(
        "text",
        x        = 4000,
        y        = 0.8,
        label    = paste("log-rank test, p.val = ", p.val),
        size     = 4.5,
        hjust    = 1,
        fontface = "bold"
      )
  )
  # rho = 1 the Gehan-Wilcoxon test
  print(DNFA)
  print(stats)
  #  fit = survfit(SurvObj ~ df$Group, data = df)
  #  tst <- comp(fit)$tests$lrTests
  #  print(tst)
}

##############################################################################
##  Kaplan-Meier curve with clinic and DNFA RNASeq data in each tumor group ##
##############################################################################
plot.km.DNFA.each.tumor <- function(DNFA, tumor) {
  DNFA.TCGA.GTEX <- read.csv(
    file.path("project-data",
      "DNFASKCMandGTEX.csv"),
    header = TRUE,
    sep = ","
  )
  DNFA.TCGA <- DNFA.TCGA.GTEX[DNFA.TCGA.GTEX$study == 'TCGA', ]
  DNFA.TCGA <-
    DNFA.TCGA[DNFA.TCGA$primary.disease.or.tissue == tumor, ]
  DNFA.TCGA <-
    DNFA.TCGA[DNFA.TCGA$sample_type != "Solid Tissue Normal",]
  #  DNFA.TCGA <- DNFA.TCGA[which(DNFA.TCGA$OS.time < 4001), ]
  #  DNFA.TCGA <- subset(DNFA.TCGA,  OS.time < 4000)
  #  DNFA.TCGA <- subset(DNFA.TCGA,  OS.time > 0)
  DNFA.TCGA <- droplevels(DNFA.TCGA)
  df <- na.omit(DNFA.TCGA)
  number <- nrow(df)
  sub <- round(number / 5, digits = 0)
  bottom.label <- paste("Bottom 20%, n = ", sub)
  top.label <- paste("Top 20%, n = ", sub)
  df$Group[df[[DNFA]] < quantile(df[[DNFA]], prob = 0.2)] = "Bottom 20%"
  df$Group[df[[DNFA]] > quantile(df[[DNFA]], prob = 0.8)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS.time, OS == 1))
  df <- na.omit(df)
  km <-
    survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <-
    survdiff(SurvObj ~ df$Group, data = df, rho = 0) # rho = 0 log-rank
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
  #  p.val.2 <- surv_pvalue(km)
  black.bold.12pt <- element_text(
    face   = "bold",
    size   = 12,
    family = "Tahoma",
    colour = "black"
  )
  print(
    ggplot2::autoplot(
      km,
      xlab = "Days",
      xlim = c(0, 4000),
      ylab = "Survival Probability",
      main = paste0("Kaplan-Meier plot of ", tumor, " (",
        number, " cases)")
    ) +
      theme_bw() +
      theme(
        plot.title           = black.bold.12pt,
        axis.title           = black.bold.12pt,
        axis.text            = black.bold.12pt,
        axis.line.x          = element_line(color  = "black"),
        axis.line.y          = element_line(color  = "black"),
        panel.grid           = element_blank(),
        strip.text           = black.bold.12pt,
        legend.text          = black.bold.12pt ,
        legend.title         = black.bold.12pt ,
        legend.position      = c(1, 1),
        legend.justification = c(1, 1)
      ) +
      guides(fill = FALSE) +
      scale_color_manual(
        values = c("red", "blue"),
        name   = paste(DNFA, "mRNA expression"),
        breaks = c("Bottom 20%", "Top 20%"),
        labels = c(bottom.label, top.label)
      ) +
      geom_point(size = 0.25) +
      annotate(
        "text",
        x        = 4000,
        y        = 0.8,
        label    = paste("log-rank test, p.val = ", p.val),
        size     = 4.5,
        hjust    = 1,
        fontface = "bold"
      )
  )
  # rho = 1 the Gehan-Wilcoxon test
  print(DNFA)
  print(stats)
  #  fit = survfit(SurvObj ~ df$Group, data = df)
  #  tst <- comp(fit)$tests$lrTests
  #  print(tst)
}

####################################################
## Boxplots of relative levels of DNFA expression ##
####################################################
plot.DNFA.TCGA.GTEX (get.DNFA.TCGA.GTEX.RNAseq.long())

plot.DNFA.TCGA (get.DNFA.TCGA.RNAseq.long())

plot.DNFA.GTEX (get.DNFA.GTEX.RNAseq.long())

plot.DNFA.TCGA.each.tumor (x = get.DNFA.TCGA.RNAseq.long(),
  y = "Prostate Adenocarcinoma")

lapply(get.disease.list(),
  plot.DNFA.TCGA.each.tumor,
  x = get.DNFA.TCGA.RNAseq.long())

###########################################################################
## PCA with DNFA or DNCS gene expression on SKCM and normal skin tissues ##
###########################################################################

DNFA.list <- c("ACLY", "ACSS2", "ACACA", "SCD", "FASN", "SREBF1")
DNFA.list1 <- c(DNFA.list, "ACSL1")
DNCS.list <- c("HMGCS1", "HMGCR", "MVK", "PMVK", "MVD", "SREBF2")
three.list <- list(DNFA.list, DNFA.list1, DNCS.list)

lapply(three.list, plot.DNFA.PCA)

#############################################################################
##  Kaplan-Meier curve with clinic and DNFA RNASeq data in all tumor group ##
#############################################################################
plot.km.DNFA.all.tumors("SCD")
DNFA.gene <- c("ACLY",
  "ACSS2",
  "ACACA",
  "SCD",
  "FASN",
  "ACSL1",
  "HMGCS1",
  "HMGCR",
  "MVK",
  "PMVK")
names(DNFA.gene) <- DNFA.gene
lapply(DNFA.gene, plot.km.DNFA.all.tumors)

##############################################################################
##  Kaplan-Meier curve with clinic and DNFA RNASeq data in each tumor group ##
##############################################################################
lapply(get.disease.list(),
  plot.km.DNFA.each.tumor,
  DNFA = "SCD")

lapply(DNFA.gene,
  plot.km.DNFA.each.tumor,
  tumor = "Bladder Urothelial Carcinoma")

lapply(DNFA.gene,
  plot.km.DNFA.each.tumor,
  tumor = "Sarcoma")

lapply(DNFA.gene,
  plot.km.DNFA.each.tumor,
  tumor = "Skin Cutaneous Melanoma")
