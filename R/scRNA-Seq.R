library(R.utils)  # needed for fread() to support.gz files directly
library(data.table)
library(ggplot2)
library(gplots)
library(grid)
library(plyr)

# TODO(dlroxe): unify this with identical function in data prep .R file.
# Return the value of the DNFA_generatedDataRoot environment variable.  If
# that variable isn't set, return "/usr/local/DNFA-genfiles/data".
.getDataDir2 <- function() {
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
.get_rnaseq_data_filename <- function() {
  local_file = file.path(
    .getDataDir2(),
    "r-extdata",
    "GSE72056_melanoma_single_cell_revised_v2.txt.gz")

    if (file.exists(local_file)) {
    return(local_file)
  }

  # TODO(dlroxe): fread() supposedly supports URLs, but this seems not to work;
  # for now it's best if the file exists locally (because it was fetched by the
  # Makefile).
  return("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72056/suppl/GSE72056_melanoma_single_cell_revised_v2.txt.gz")
}

# Select specific lipogenesis genes from data
# This function modifies it as follows:
#
# 1. Transpose the data so that cells are rows rather than columns.
# 2. Create a data frame from transposed data.
# 3. Rename a column that describes malignant state to 'malignancy'.
# 4. Renames numeric elements of the "malignancy" column with the
#    English words "Unresolved", "Malignant", or "Non-malignant".
#
# The resulting data frame is returned to the caller.
#
get_lipogenesis_data <- function(data) {
  # -----------------------------------------------------------------------------------------------
  #select the genes of interest
  SREBF1 <- data[data$Cell == "SREBF1",]
  SREBF2 <- data[data$Cell == "SREBF2",]
  FASN <- data[data$Cell == "FASN",]
  SCD <- data[data$Cell == "SCD",]
  ACACA <- data[data$Cell == "ACACA",]
  ACSS2 <- data[data$Cell == "ACSS2",]
  ACSL1 <- data[data$Cell == "ACSL1",]
  ACLY <- data[data$Cell == "ACLY",]
  HMGCR <- data[data$Cell == "HMGCR",]
  HMGCS1 <- data[data$Cell == "HMGCS1",]
  PPARGC1A <- data[data$Cell == "PPARGC1A",]
  PPARGC1B <- data[data$Cell == "PPARGC1B",]
  MITF <- data[data$Cell == "MITF",]
  AXL <- data[data$Cell == "AXL",]
  tumor <- data[data$Cell == "tumor",]
  malignant <- data[data$Cell == "malignant(1=no,2=yes,0=unresolved)",]

  # -----------------------------------------------------------------------------------------------
  #combine them into a new table "geneset"
  geneset <- rbind(tumor,malignant,
                   SREBF1,SREBF2,FASN,SCD,ACACA,ACSS2,ACLY,ACSL1,
                   HMGCR,HMGCS1,PPARGC1A,PPARGC1B,MITF,AXL)

  ## transpose the table
  # first remember the names
  n <- geneset$Cell
  # transpose all but the first column (name)
  transposed_geneset <- as.data.frame(t(geneset[,-1]))
  colnames(transposed_geneset) <- n
  renamed_geneset <- rename(
    transposed_geneset,
    c('malignant(1=no,2=yes,0=unresolved)' = 'malignancy'))
  renamed_geneset$malignancy <- factor(renamed_geneset$malignancy)
  renamed_geneset$malignancy <- revalue(renamed_geneset$malignancy,
                                        c("0" = "Unresolved",
                                          "1" = "Non-malignant",
                                          "2" = "Malignant"))
  return(renamed_geneset)
}


# Single cell RNA-seq data of melanomas were downloaded from GSE72056.
# the downloaded txt file GSE72056_melanoma_single_cell_revised_v2.txt is
# produced locally by the Makefile.
.get_rnaseq_data <- function() {
  # The data table is over 300 Mb, so it seems better to use the fread
  # function in R package data.table to quickly import the dataset as a dataframe.
  melanomaSingleCellFile <- .get_rnaseq_data_filename()
  singleRNAseq <- fread(melanomaSingleCellFile, header = T)
  sampleSingleRNAseq <- head(singleRNAseq)[,c(1,2,3,4,5,6,7)]
  View(sampleSingleRNAseq)
  return(singleRNAseq)
}

# select RNA-seq data from malignant or non-malignant cells
.get_totalgeneset <- function(lipogenesis_data) {
  totalgeneset <- subset(lipogenesis_data, malignancy != "Unresolved")
  totalgeneset$tumor <- as.factor(totalgeneset$tumor)
  levels(totalgeneset$tumor)
  return(totalgeneset)
}

make_single_cell_plots <- function() {

  singleRNAseq <- .get_rnaseq_data()
  lipogenesis_data <- get_lipogenesis_data(singleRNAseq)
  totalgeneset <- .get_totalgeneset(lipogenesis_data)
  # malignantgeneset <- subset(lipogenesis_data, malignancy == "Malignant")
  # nonmalignantgeneset <- subset(lipogenesis_data, malignancy == "Non-malignant")

  singleCellBoxplot <- geom_boxplot(
    size = 1, outlier.colour = NA, color = "black", aes(fill = malignancy))

  singleCellTheme <-
    theme_bw() +
    theme(
      axis.title = element_text(face = "bold", size = 18, color = "black"),
      axis.text = element_text(size = 18,face = "bold", color = "black"),
      axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
      axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
      axis.line.x = element_line(color = "black", size = 1),
      axis.line.y = element_line(color = "black", size = 1),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(.25, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.text = element_text(size = 18, face = "bold", colour = 'black'),
      legend.position = c(0,1),
      legend.justification = c(-0.1,1.1),
      legend.key.height = unit(2.2, 'lines'),
      legend.background = element_rect(fill = "white"))

  singleCellGuidesAndScales <-
    guides(fill = guide_legend(title = NULL)) +
    scale_y_continuous(breaks = seq(0, 10, 2)) +
    scale_fill_discrete(labels = c("Nonmalignant cells", "Malignant cells"))

  ###############################################################################
  # compare single cell SREBF1 expression level in malignant and non-malignant
  # cells using boxplot graph within the same graph##
  ###############################################################################
  plot1 <- ggplot(totalgeneset, aes(x = factor(tumor), y = SREBF1)) +
    singleCellBoxplot + singleCellTheme + singleCellGuidesAndScales +
    labs(x = "tumor samples", y = "SREBF1 mRNA counts")

  print(plot1)
  # -----------------------------------------------------------------------------------------------
  plot2 <- ggplot(totalgeneset, aes(x = factor(tumor), y = FASN)) +
    singleCellBoxplot + singleCellTheme + singleCellGuidesAndScales +
    labs(x = "tumor samples", y = "FASN mRNA counts")

  print(plot2)
  # -----------------------------------------------------------------------------------------------

  plot3 <- ggplot(totalgeneset, aes(x=factor(tumor), y = SCD)) +
    singleCellBoxplot + singleCellTheme + singleCellGuidesAndScales +
    labs(x = "tumor samples", y = "SCD mRNA counts")

  print(plot3)
  # -----------------------------------------------------------------------------------------------
  plot4 <- ggplot(totalgeneset, aes(x=factor(tumor), y = ACACA)) +
    singleCellBoxplot + singleCellTheme + singleCellGuidesAndScales +
    labs(x = "tumor samples", y = "ACACA mRNA counts")

  print(plot4)
  # -----------------------------------------------------------------------------------------------
  plot5 <- ggplot(totalgeneset, aes(x=factor(tumor), y = SREBF2)) +
    singleCellBoxplot + singleCellTheme + singleCellGuidesAndScales +
    labs(x = "tumor samples", y = "SREBF2 mRNA counts")

  print(plot5)
  # -----------------------------------------------------------------------------------------------
  plot6 <- ggplot(totalgeneset, aes(x=factor(tumor), y = MITF)) +
    singleCellBoxplot + singleCellTheme + singleCellGuidesAndScales +
    labs(x = "tumor samples", y = "MITF mRNA counts")

  print(plot6)
  # -----------------------------------------------------------------------------------------------
  plot7 <- ggplot(totalgeneset, aes(x=factor(tumor), y = AXL)) +
    singleCellBoxplot + singleCellTheme + singleCellGuidesAndScales +
    labs(x = "tumor samples", y = "AXL mRNA counts")

  print(plot7)
}

make_single_cell_plots()
