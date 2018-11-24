# The following R script combines all gene count results of each samples from RNA-Seq experiment and generate a mega data set for DESEQ analysis in the seperate R scripts. 
# Combine a raw data file per sample to a single merged file.

# Sample format of each file "testN-KReadsPerGene.out.tab":
# N_unmapped      105984  105984 105984
# N_multimapping  171224  171224 171224
# N_noFeature     177533 3433150 277677
# N_ambiguous     319796    9239 136511
# ENSG00000223972      0       0      0
# ENSG00000227232     10       0     10
# ENSG00000278267      0       0      0
# ...
# STAR outputs read counts per gene into ReadsPerGene.out.tab with 4 columns
# column 1: gene ID
# column 2: counts for unstranded RNA-seq
# column 3: counts for the 1st read strand aligned with RNA
# column 4: counts for the 2nd read strand aligned with RNA
#
# Our RNA-Seq libraries were contructed by the
# NEBNext® Ultra™ Directional RNA Library Prep Kit for Illumina.
# This kit uses dUTP method to generate anti-sense strand for the 1st read strand synthesis.
# (the original RNA strand is degradated due to the dUTP incorporated),
# so the 2nf read strand in the result column is from the original RNA strand.
# See the illustration from the following link https://bit.ly/29Yi771
# We choose the data from columns 4 as the gene counts for the given sample.
#
#
# The notional desired format of the merged data is:
#                 test1.1.1 test1.1.2 test1.1.3 test1.2.1 test1.2.2 ...
# N_unmapped      105984       105984    105984
# N_multimapping  171224       171224    171224
# N_noFeature     177533      3433150    277677
# N_ambiguous     319796         9239    136511
# ENSG00000223972      0            0        0
# ENSG00000227232     10            0       10
# ENSG00000278267      0            0        0
# ...
#
# However, only the .3 columns are meaningful for the
# analyses to be performed by this package.  Therefore,
# the final merged output should look like this:
#
#                 test1.1.3 test1.2.3 test1.3.3 test1.4.3 test2.1.3 ...
# N_unmapped
# N_multimapping
# N_noFeature
# N_ambiguous
# ENSG00000223972
# ENSG00000227232
# ENSG00000278267
# ...
#

#
# Row names indicate metadata like the number of unmapped segments,
# and Ensembl gene designations.  The designations are translated
# to human-readable form and appended as a final column in the merged
# data, taking the first symbol (if any) for each designator, using
# the "org.Hs.eg.db" package:
#                 ...    SYMBOL
# N_unmapped      ...      null
# N_multimapping  ...      null
# N_noFeature     ...      null
# N_ambiguous     ...      null
# ENSG00000223972 ...   DDX11L1
# ENSG00000227232 ...    WASH7P
# ENSG00000278267 ... MIR6859-1
# ...

library("AnnotationDbi")
library("base")
library("dplyr")
library("org.Hs.eg.db")
library("utils")

makeColumnBaseNameFromFileName <- function(fileName) {
  # Generate 'test6.1.' from 'test6_S4_L001ReadsPerGene.out.tab'.
  newColBaseName <- paste(substring(fileName, 1, 5),
                          substring(fileName, 13, 13),
                          '',  # gives us a trailing '.'
                          sep = '.')
  print(paste("name", newColBaseName))
  return(newColBaseName)
}

# Reframe a data file with desired row and column names.
processReadsFile <- function(fullFileName) {
  print(paste("processReadsFile(): ", fullFileName))
  table <- read.table(fullFileName, stringsAsFactors = T)

  # Reframe the data, taking row names from first column.
  table.with.rownames <- data.frame(table[,-1], row.names = table[,1])

  newColBaseName <- makeColumnBaseNameFromFileName(basename(fullFileName))
  
  # Rename the columns as test6.4.1, test6.4.2, test6.4.3.
  colnames(table.with.rownames) <- paste0(newColBaseName, 1:3)

  ## ...but then preserve only column 3
  name3 <- paste(newColBaseName, 3, sep = '')
  table.with.rownames <- dplyr::select(table.with.rownames, name3)
  print(head(table.with.rownames))
  return(table.with.rownames)
}

# Merge reframed data from all input files, and add 'symbol' column.
mergeFiles <- function(files) {
  processedData <- lapply(files, processReadsFile)
  test <- do.call(cbind.data.frame, processedData)

  test$symbol <- mapIds(org.Hs.eg.db,
                        keys = row.names(test),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")
  return(test)
}

.getDataDir <- function() {
  return(Sys.getenv("DNFA_generatedDataRoot", unset = "/usr/local/DNFA-genfiles/data"))
}

.getFile <- function(stem) {
  name <- paste('test', stem, 'ReadsPerGene.out.tab', sep = "")
  return(file.path(.getDataDir(), name))
}

prepare_data <- function() {
  stems <- c(
    "1_S2_L001", "1_S2_L002", "1_S2_L003", "1_S2_L004",
    "2_S3_L001", "2_S3_L002", "2_S3_L003", "2_S3_L004",
    "3_S4_L001", "3_S4_L002", "3_S4_L003", "3_S4_L004",
    "4_S6_L001", "4_S6_L002", "4_S6_L003", "4_S6_L004",
    "5_S5_L001", "5_S5_L002", "5_S5_L003", "5_S5_L004",
    "6_S1_L001", "6_S1_L002", "6_S1_L003", "6_S1_L004"
  )

  write.csv(mergeFiles(lapply(stems, .getFile)),
            file = file.path(.getDataDir(), "r-extdata", 'gene-counts.csv'))
}
