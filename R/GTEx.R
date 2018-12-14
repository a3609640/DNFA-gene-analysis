library(data.table)
library(ggplot2)

# TODO(dlroxe): unify this with identical function in data prep .R file.
# Return the value of the DNFA_generatedDataRoot environment variable.  If
# that variable isn't set, return "/usr/local/DNFA-genfiles/data".
.getDataDir3 <- function() {
  return(Sys.getenv("DNFA_generatedDataRoot", unset = "/usr/local/DNFA-genfiles/data"))
}

# In case of difficulty coordinating with Makefile output, this function
# can be adjusted to return some other string, e.g. return("~/foo.txt.gz").
# In theory, we could return a URL to a .gz file here as well, but then it
# would be difficult to ensure that the file remains available locally for
# repeated invocations.  It should be possible to execute this program
# many times after downloading the file only once.  A possible alternative
# is to check the existence of the local file, then resort to a URL
# if the local file is missing.
.get_gtex_data_filename <- function() {
  local_file = file.path(
    .getDataDir3(),
    "r-extdata",
    "GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz")

  if (file.exists(local_file)) {
    return(local_file)
  }

  # TODO(dlroxe): fread() supposedly supports URLs, but this seems not to work;
  # for now it's best if the file exists locally (because it was fetched by the
  # Makefile).
  return("http://storage.googleapis.com/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz")
}

.get_gtex_annotations_filename <- function() {
  local_file = file.path(
    .getDataDir3(),
    "r-extdata",
    "GTEx_Data_V6_Annotations_SampleAttributesDS.txt")
  if (file.exists(local_file)) {
    return(local_file)
  }

  # TODO(dlroxe): fread() supposedly supports URLs, but this seems not to work;
  # for now it's best if the file exists locally (because it was fetched by the
  # Makefile).
  return("http://storage.googleapis.com/gtex_analysis_v6p/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt")
}

# Looks like this site might be a useful pointer towards a source for the data:
# https://github.com/joed3/GTExV6PRareVariation
## showProgress = T is necessary, 
## otherwise "Error: isLOGICAL(showProgress) is not TRUE"
.get_gene_rpkm <- function() {
  gene <- data.table::fread(.get_gtex_data_filename(), header = T, showProgress = T)
  return(gene)
}

.get_annotations <- function() {
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

doGTEX <- function() {
  gene_rpkm <- .get_gene_rpkm()
  Annotations <- .get_annotations()
  plot1 <- .plot_goi(gene_rpkm, Annotations, goi = "FASN")
  plot2 <- .plot_goi(gene_rpkm, Annotations, goi = "SCD")
  plot3 <- .plot_goi(gene_rpkm, Annotations, goi = "SREBF1")
  plot4 <- .plot_goi(gene_rpkm, Annotations, goi = "HMGCR")
  plot5 <- .plot_goi(gene_rpkm, Annotations, goi = "HMGCS1")
  plot6 <- .plot_goi(gene_rpkm, Annotations, goi = "SREBF2")
  plot7 <- .plot_goi(gene_rpkm, Annotations, goi = "MITF")
  print(plot1)
  print(plot2)
  print(plot3)
  print(plot4)
  print(plot5)
  print(plot6)
  print(plot7)
}

doGTEX()
