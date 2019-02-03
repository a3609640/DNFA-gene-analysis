library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)

# TODO(dlroxe): unify this with identical function in data prep .R file.
# Return the value of the DNFA_generatedDataRoot environment variable.  If
# that variable isn't set, return "/usr/local/DNFA-genfiles/data".
.getDataDir3 <- function() {
  return(Sys.getenv(
    "DNFA_generatedDataRoot",
    unset = file.path(Sys.getenv("HOME"), "Downloads")
  ))
  }

# Returns the fully-qualified local filename for GTEx data.  If the
# file cannot be found at an expected location, fetches the data and
# stores it there.  (In theory the Makefile should download and store
# the file before this script is run; but perhaps it makes more sense
# to take that responsibility away from the Makefil and so avoid any
# difficulty coordinating between its output and this script's input.)
.get_gtex_data_filename <- function() {
  local_file = file.path(
    .getDataDir3(),
    # "r-extdata",  # TODO(dlroxe): probably, stop using r-extdata in Makefile
#    "GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz")
    "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz")
  
  if (!file.exists(local_file)) {
    download.file(
#      url = "http://storage.googleapis.com/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz",
      url = "https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz",
      destfile = local_file)
  }
  
  return(local_file)
}

.get_gtex_annotations_filename <- function() {
  local_file = file.path(
    .getDataDir3(),
    # "r-extdata",  # TODO(dlroxe): probably, stop using r-extdata in Makefile
 #   "GTEx_Data_V6_Annotations_SampleAttributesDS.txt")
    "GTEx_v7_Annotations_SampleAttributesDS.txt")

  if (!file.exists(local_file)) {
    download.file(
 #     url = "http://storage.googleapis.com/gtex_analysis_v6p/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt",
      url = "https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt",
      destfile = local_file)
  }
  
  return(local_file)
}

# Looks like this site might be a useful pointer towards a source for the data:
# https://github.com/joed3/GTExV6PRareVariation
## showProgress = T is necessary, 
## otherwise "Error: isLOGICAL(showProgress) is not TRUE"
.get_gene <- function() {
  gene <- data.table::fread(.get_gtex_data_filename(), header = T, showProgress = T)
  return(gene)
}

.get_gene_annotations <- function() {
  Annotations <- data.table::fread(
    .get_gtex_annotations_filename(), header = T, showProgress = T)
  Annotations <- Annotations[,c(1,7)]
  return(Annotations)
}


.plot_goi <- function(goi, gene, gene_annotations) {
  go <- gene[gene$Description == goi,]
  go <- go[, -c(1,2)]
  go <- t(go)
  ## go is generated as a matrix, and it has to be converted into data frame.
  go <- data.frame(go)
  setDT(go, keep.rownames = TRUE)[]
  colnames(go) <- c("SAMPID", "goi")
  ## one line option is: df$names<-rownames(df)
  goexpression <- merge(go, gene_annotations, by = 'SAMPID')
  ## somehow the numbers of SREBF1 columns are all changed into character
  goexpression$SMTSD <- as.factor(goexpression$SMTSD)
  #  In particular, as.numeric applied to a factor is meaningless,
  #  and may happen by implicit coercion.
  #  To transform a factor f to approximately its original numeric values,
  #  as.numeric(levels(f))[f] is recommended.
#  goexpression$goi <- as.numeric(levels(goexpression$goi))[goexpression$goi]
  goexpression <- as.data.frame(goexpression)
  ## draw boxplot for FASN expression across different tissues
  mean <- within(goexpression, SMTSD <- reorder(SMTSD, log10(goi), median))
  black.bold.12pt <- ggplot2::element_text(face   = "bold",
                                           size   = 12,
                                           colour = "black")
  genePlot <- ggplot(mean, aes(x = SMTSD, y = log10(goi))) +
   geom_boxplot() + theme_bw() +
   labs(x = "Tissue types (GTEx)",
        y = paste("log10(", goi, "RNA counts)")) +
        theme(axis.title           = black.bold.12pt,
              axis.text            = ggplot2::element_text(size  = 12,
                                                           angle = 90,
                                                           hjust = 1,
                                                           face  = "bold",
                                                           color = "black"),
              axis.line.x          = ggplot2::element_line(color = "black"),
              axis.line.y          = ggplot2::element_line(color = "black"),
              panel.grid           = ggplot2::element_blank(),
              strip.text           = black.bold.12pt,
              legend.text          = black.bold.12pt,
              legend.title         = black.bold.12pt,
              legend.justification = c(1,1)) +
              guides(fill = guide_legend(title = NULL))
              print(genePlot)
}

##################################################################
## important to keep gene and gene_annotations as global variable. They are too big in size.
## make sure read them once and avoid to construc functions over them.



get.EIF.RNAseq.GTEx <- function(gene, gene_annotations){
  ## use %in% instead of == for subsetting rows!!
  EIF.gene <- c("EIF4A1","EIF4B","EIF4E","EIF4G1","EIF4EBP1","RPS6KB1","MYC")
  EIF <- gene[gene$Description %in% EIF.gene,]
  EIF <- EIF[, -1]
  n <- EIF$Description
  ## go is generated as a matrix, and it has to be converted into data frame.
  EIF <- as.data.frame(t(EIF[,-1]))
  colnames(EIF) <- n
  setDT(EIF, keep.rownames = TRUE)[]
# colnames(EIF) <- c("SAMPID", "goi")
  colnames(EIF) <- c("SAMPID", n)
## one line option is: df$names<-rownames(df)
  EIF.RNAseq.GTEx <- merge(EIF, gene_annotations, by = 'SAMPID')
## somehow the numbers of SREBF1 columns are all changed into character
  EIF.RNAseq.GTEx$SMTSD <- as.factor(EIF.RNAseq.GTEx$SMTSD)
  sapply(EIF.RNAseq.GTEx, class)
  EIF.RNAseq.GTEx <- as.data.frame(EIF.RNAseq.GTEx[,-1])
  tissues <- levels(EIF.RNAseq.GTEx$SMTSD)
  class(tissues)
  return(EIF.RNAseq.GTEx)
  }

get.EIF.score.GTEx <- function(gene, gene_annotations){
# rm(EIFscore)
  EIF.RNAseq.GTEx <- get.EIF.RNAseq.GTEx(gene, gene_annotations)
  EIF.score.GTEx <- EIF.RNAseq.GTEx
  EIF.score.GTEx$EIF4A1score <- EIF.RNAseq.GTEx$EIF4A1/EIF.RNAseq.GTEx$EIF4E
  EIF.score.GTEx$EIF4Escore <- EIF.RNAseq.GTEx$EIF4E/EIF.RNAseq.GTEx$EIF4E
  EIF.score.GTEx$EIF4G1score <- EIF.RNAseq.GTEx$EIF4G1/EIF.RNAseq.GTEx$EIF4E
  EIF.score.GTEx$EIF4EBP1score <- EIF.RNAseq.GTEx$EIF4EBP1/EIF.RNAseq.GTEx$EIF4E
  EIF.score.GTEx$RPS6KB1score <- EIF.RNAseq.GTEx$RPS6KB1/EIF.RNAseq.GTEx$EIF4E
#  EIF.score.GTEx <- EIF.score.GTEx [, c("EIF4A1score","EIF4Escore",
#                                        "EIF4G1score","EIF4EBP1score",
#                                        "RPS6KB1score")]
  return(EIF.score.GTEx)
  }


############################################
###  plot RNAseq from all tissue samples ###
############################################

plotEIF <-  function (x) {
  name <- deparse(substitute(x))
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

plot.EIFandScore.all.tissues <- function (gene, gene_annotations){
  EIF.RNAseq.GTEx <- get.EIF.RNAseq.GTEx(gene, gene_annotations)
  EIF.RNAseq.GTEx <- EIF.RNAseq.GTEx[, 
                                     c("EIF4A1","EIF4E","EIF4G1", 
                                       "EIF4EBP1","RPS6KB1","SMTSD")]
  EIF.score.GTEx <- get.EIF.score.GTEx(gene, gene_annotations)
  EIF.score.GTEx <- EIF.score.GTEx[, c("EIF4Escore","EIF4A1score",
                                       "EIF4G1score","EIF4EBP1score",
                                       "RPS6KB1score")]
  number <- nrow(EIF.RNAseq.GTEx)
  EIF.RNAseq.GTEx.all.tissues <- melt(EIF.RNAseq.GTEx)
  EIF.score.GTEx.all.tissues <- melt(EIF.score.GTEx)
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
  p1 <- plotEIF(EIF.RNAseq.GTEx.all.tissues) +
    labs(title = paste0("All healthy tissues n = ", number),
         x     = "eIF4F subunit RNAseq",
         y     = paste0("log2(value)")) +
    stat_compare_means(comparisons = my_comparison1, method = "t.test")
  p1$layers[[2]]$aes_params$textsize <- 5  
  p2 <- plotEIF(EIF.score.GTEx.all.tissues) + 
    labs(title = paste0("All healthy tissues n = ", number),
         x     = "eIF4E ratio score",
         y     = paste0("log2(value)")) +
    stat_compare_means(comparisons = my_comparison2, method = "t.test")
  p2$layers[[2]]$aes_params$textsize <- 5
  grid.arrange(p1, p2, ncol = 2)
  }



plot.EIFandScore.each.tissue <-
  function(m, gene, gene_annotations) {
    gene_names <- c("EIF4A1", "EIF4E", "EIF4G1", "EIF4EBP1", "RPS6KB1", "SMTSD")
    gene_score_names <- c(
      "EIF4A1score",
      "EIF4Escore",
      "EIF4G1score",
      "EIF4EBP1score",
      "RPS6KB1score",
      "SMTSD"
    )
  EIF.RNAseq.GTEx <- get.EIF.RNAseq.GTEx(gene, gene_annotations)
  EIF.RNAseq.GTEx <- EIF.RNAseq.GTEx[, gene_names]
  EIF.RNAseq.GTEx <- EIF.RNAseq.GTEx[EIF.RNAseq.GTEx$SMTSD == m,]

  EIF.score.GTEx <- get.EIF.score.GTEx(gene, gene_annotations)
  EIF.score.GTEx <- EIF.score.GTEx[, gene_score_names]
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

  if (medianEIF4G1score >= medianEIF4EBP1score) {
    print(paste("EIF is activated in", m))
  } else {
    print(paste("EIF is inhibited in", m))
  }

  par(mfrow = c(1, 2))
  p1 <- plotEIF(EIF.RNAseq.GTEx.each.tissues) +
    labs(
      title = paste0(m, " n = ", number),
      x     = "eIF4F subunit RNAseq",
      y     = paste0("log2(value)")
    ) +
    stat_compare_means(comparisons = my_comparison1, method = "t.test")
  p1$layers[[2]]$aes_params$textsize <- 5
  
  p2 <- plotEIF(EIF.score.GTEx.each.tissues) +
    labs(
      title = paste0(m, " n = ", number),
      x     = "eIF4E ratio score",
      y     = paste0("log2(value)")
    ) +
    stat_compare_means(comparisons = my_comparison2, method = "t.test")
  p2$layers[[2]]$aes_params$textsize <- 5
  grid.arrange(p1, p2, ncol = 2)
}


gene <- .get_gene()
gene_annotations <- .get_gene_annotations()

.plot_goi("EIF4EBP1", gene, gene_annotations)

EIF.gene <- c("EIF4A1","EIF4E","EIF4G1","EIF4EBP1","RPS6KB1","MYC")
names(EIF.gene) <- EIF.gene

DNFA.gene <- c("SCD", "FASN", "ACLY","ACSS2","SREBF1",
               "HMGCR","HMGCS1","SREBF2","MITF")
names(DNFA.gene ) <- DNFA.gene
lapply(EIF.gene, .plot_goi, gene=gene, gene_annotations=gene_annotations)

plot.EIFandScore.all.tissues(gene=gene, gene_annotations=gene_annotations)

plot.EIFandScore.each.tissue("Muscle - Skeletal", 
                             gene=gene, 
                             gene_annotations=gene_annotations)
EIF.RNAseq.GTEx <- get.EIF.RNAseq.GTEx(gene=gene, 
                                       gene_annotations=gene_annotations)

tissues <- levels(EIF.RNAseq.GTEx$SMTSD)
sapply(tissues, plot.EIFandScore.each.tissue, 
       gene=gene, gene_annotations=gene_annotations)








                        