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
    "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz")
  
  if (!file.exists(local_file)) {
    download.file(
#      url = "http://storage.googleapis.com/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz",
      url = "https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz",
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

.plot_goi <- function(gene, Annotations, goi) {
  go <- gene[gene$Description == goi,]
  go <- t(go)
  ## go is generated as a matrix, and it has to be converted into data frame.
  go <- data.frame(go)
  setDT(go, keep.rownames = TRUE)[]
  colnames(go) <- c("SAMPID", "goi")
  ## one line option is: df$names<-rownames(df)
  goexpression <- merge(go, Annotations, by = 'SAMPID')
  ## somehow the numbers of SREBF1 columns are all changed into character
  goexpression$SMTSD <- as.factor(goexpression$SMTSD)
  #  In particular, as.numeric applied to a factor is meaningless,
  #  and may happen by implicit coercion.
  #  To transform a factor f to approximately its original numeric values,
  #  as.numeric(levels(f))[f] is recommended.
  goexpression$goi <- as.numeric(levels(goexpression$goi))[goexpression$goi]
  goexpression <- as.data.frame(goexpression)
  ## draw boxplot for FASN expression across different tissues
  mean <- within(goexpression, SMTSD <- reorder(SMTSD, log2(goi), median))
  black.bold.12pt <- ggplot2::element_text(face   = "bold",
                                           size   = 12,
                                           colour = "black")
  genePlot <- ggplot(mean, aes(x = SMTSD, y = log2(goi))) +
   geom_boxplot() + theme_bw() +
   labs(x = "Tissue types (GTEx)",
        y = paste("log2(", goi, "RNA counts)")) +
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

plotGOI_DNFA <- function(gene, annotations) {
  plot1 <- .plot_goi(gene, annotations, goi = "FASN")
  plot2 <- .plot_goi(gene, annotations, goi = "SCD")
  plot3 <- .plot_goi(gene, annotations, goi = "SREBF1")
  plot4 <- .plot_goi(gene, annotations, goi = "HMGCR")
  plot5 <- .plot_goi(gene, annotations, goi = "HMGCS1")
  plot6 <- .plot_goi(gene, annotations, goi = "SREBF2")
  plot7 <- .plot_goi(gene, annotations, goi = "MITF")
  print(plot1)
  print(plot2)
  print(plot3)
  print(plot4)
  print(plot5)
  print(plot6)
  print(plot7)
}

plotGOI_EIF <- function(gene, annotations) {
  plot1 <- .plot_goi(gene, annotations, goi = "EIF4A1")
  plot2 <- .plot_goi(gene, annotations, goi = "EIF4E")
  plot3 <- .plot_goi(gene, annotations, goi = "EIF4G1")
  plot4 <- .plot_goi(gene, annotations, goi = "EIF4EBP1")
  plot5 <- .plot_goi(gene, annotations, goi = "RPS6KB1")
  plot6 <- .plot_goi(gene, annotations, goi = "MYC")
  print(plot1)
  print(plot2)
  print(plot3)
  print(plot4)
  print(plot5)
  print(plot6)
}

## TO do: divide different tissues into groups by comparing means of 
## EIF4Gscore to EIF4EBP1score





##################################################################
### TO DO: organize the following scripts , (they worked)


# plotGOI_DNFA(gene, gene_annotations)
# plotGOI_EIF(gene, gene_annotations)
## important to keep gene and gene_annotations as global variable. They are too big in size.
## make sure read them once and avoid to construc functions over them.
gene <- .get_gene()
gene_annotations <- .get_gene_annotations()

get.EIF.RNAseq.GTEx <- function(){
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
  colnames(EIF) <- c("SAMPID", EIF.gene)
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

get.EIF.score.GTEx <- function(){
# rm(EIFscore)
  EIF.RNAseq.GTEx <- get.EIF.RNAseq.GTEx()
  EIF.score.GTEx <- EIF.RNAseq.GTEx
  EIF.score.GTEx$EIF4Escore <- EIF.RNAseq.GTEx$EIF4E/EIF.RNAseq.GTEx$EIF4E
  EIF.score.GTEx$EIF4G1score <- EIF.RNAseq.GTEx$EIF4G1/EIF.RNAseq.GTEx$EIF4E
  EIF.score.GTEx$EIF4EBP1score <- EIF.RNAseq.GTEx$EIF4EBP1/EIF.RNAseq.GTEx$EIF4E
  EIF.score.GTEx$RPS6KB1score <- EIF.RNAseq.GTEx$RPS6KB1/EIF.RNAseq.GTEx$EIF4E
  EIF.score.GTEx <- EIF.score.GTEx [, 8:12]
  return(EIF.score.GTEx)
  }

plot.EIF.RNAseq.score <- function (m) {
  EIF.RNAseq.GTEx <- get.EIF.RNAseq.GTEx()
  EIF.score.GTEx <- get.EIF.score.GTEx()
  EIF.RNAseq.GTEx <- EIF.RNAseq.GTEx[EIF.RNAseq.GTEx$SMTSD == m,]
  EIF.score.GTEx <- EIF.score.GTEx[EIF.score.GTEx$SMTSD == m,]
  medianEIF4G1score <- median(EIF.score.GTEx$EIF4G1score)
  medianEIF4EBP1score <- median(EIF.score.GTEx$EIF4EBP1score)
  # tissue$SMTSD <- NULL
  if (medianEIF4G1score < medianEIF4EBP1score ) {
    par(mfrow=c(1,2))
    boxplot(log2(EIF.RNAseq.GTEx[, 
                               c("EIF4E", "EIF4G1", "EIF4EBP1", "RPS6KB1")]), 
            main= paste0("EIF RNAseq counts in ",m),
            las = 2)
    boxplot(log2(EIF.score.GTEx[,
                          c("EIF4Escore","EIF4G1score","EIF4EBP1score","RPS6KB1score")]),
            main= paste0("EIF scores in ", m),
            las = 2)
    print(paste("EIF is inhibited in", m))
  } else {
    print(paste("EIF is activated in", m))
  }
}

plot.EIF.RNAseq.score ("Muscle - Skeletal")
EIF.RNAseq.GTEx <- get.EIF.RNAseq.GTEx()
tissues <- levels(EIF.RNAseq.GTEx$SMTSD)
sapply(tissues, plot.EIF.RNAseq.score)

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

plot.EIFandScore.all.tissues <- function (){
  EIF.RNAseq.GTEx <- get.EIF.RNAseq.GTEx()
  EIF.RNAseq.GTEx <- EIF.RNAseq.GTEx[, 
                                     c("EIF4E", "EIF4G1", "EIF4EBP1", "RPS6KB1", "SMTSD")]
  EIF.score.GTEx <- get.EIF.score.GTEx()
  number <- nrow(EIF.RNAseq.GTEx)
  EIF.RNAseq.GTEx.all.tissues <- melt(EIF.RNAseq.GTEx)
  EIF.score.GTEx.all.tissues <- melt(EIF.score.GTEx)
  my_comparison1 <- list( c("EIF4E", "EIF4G1"), 
                          c("EIF4G1", "EIF4EBP1"), 
                          c("EIF4E", "EIF4EBP1"),
                          c("EIF4E", "RPS6KB1"))
  my_comparison2 <- list( c("EIF4Escore", "EIF4G1score"), 
                          c("EIF4G1score", "EIF4EBP1score"), 
                          c("EIF4Escore", "EIF4EBP1score"),
                          c("EIF4Escore", "RPS6KB1score"))
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
  grid.arrange(p1, p2, ncol=2)
  }
plot.EIFandScore.all.tissues()


plot.EIFandScore.each.tissues <- function (m){
  EIF.RNAseq.GTEx <- get.EIF.RNAseq.GTEx()
  EIF.RNAseq.GTEx <- EIF.RNAseq.GTEx[, 
                                     c("EIF4E", "EIF4G1", "EIF4EBP1", "RPS6KB1", "SMTSD")]
  EIF.score.GTEx <- get.EIF.score.GTEx()
  EIF.RNAseq.GTEx <- EIF.RNAseq.GTEx[EIF.RNAseq.GTEx$SMTSD == m,]
  EIF.score.GTEx <- EIF.score.GTEx[EIF.score.GTEx$SMTSD == m,]
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
}

sapply(tissues, plot.EIFandScore.each.tissues)
plot.EIFandScore.each.tissues("Muscle - Skeletal")








plotEIF(EIF.RNAseq.GTEx.all.tissues) + 
  stat_compare_means(method = "anova", label.y = 12) + # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "EIF4E") # Pairwise comparison against EIF4E
plotEIF(EIFScore)  
grid.arrange(plotEIF(RNAcounts), plotEIF(EIFScore), ncol=2)




                        