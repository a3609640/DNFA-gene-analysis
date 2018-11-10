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
# 
#
# The desired format of the merged data is:
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

library("org.Hs.eg.db")

# Reframe a data file with desired row and column names.
processFile <- function(fullFileName) {
  fileName <- basename(fullFileName)
  table <- read.table(fullFileName, stringsAsFactors=T)
  
  # Reframe the data, taking row names from first column.
  table.with.rownames <- data.frame(table[,-1], row.names=table[,1])

  # Generate 'test6.4.' from 'test6_S4_L001ReadsPerGene.out.tab'.
  newColBaseName <- paste(substring(fileName, 1, 5),
                          substring(fileName, 8, 8),
                          '',  # gives us a trailing '.'
                          sep='.')
  
  # Rename the columns as test6.4.1, test6.4.2, test6.4.3.
  colnames(table.with.rownames) <- paste0(newColBaseName, 1:3)
  return(table.with.rownames)
}

# Merge reframed data from all input files, and add 'symbol' column.
mergeFiles <- function(files) {
  processedData <- lapply(inputFiles, processFile)
  test <- do.call(cbind.data.frame, processedData)
  
  test$symbol <- mapIds(org.Hs.eg.db,
                        keys=row.names(test),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
  return(test)
}

stems <- c(
  "1_S2_L001", "1_S2_L002", "1_S2_L003", "1_S2_L004",
  "2_S3_L001", "2_S3_L002", "2_S3_L003", "2_S3_L004",
  "3_S4_L001", "3_S4_L002", "3_S4_L003", "3_S4_L004",
  "4_S6_L001", "4_S6_L002", "4_S6_L003", "4_S6_L004",
  "5_S5_L001", "5_S5_L002", "5_S5_L003", "5_S5_L004",
  "6_S1_L001", "6_S1_L002", "6_S1_L003", "6_S1_L004"
)

# TODO(dlroxe): Attempts to use 'extdata' are experimental for the moment.
#getSystemFile <- function(name) {
#  system.file("extdata", name, package = "DNFAGeneAnalysis", mustWork = TRUE)
#}
#inputFiles <- lapply(names, getSystemFile)

# TODO(dlroxe): Get rid of dataDir.  The script should just work without
# having to learn the environment configuration that was provided to the
# Makefile.
# 
# One way to do this is to have 'make' deposit data into a well-known
# location like 'inst/extdata'.  However R studio doesn't seem happy
# with the amount of data it must manage in that case (at least, for
# all generated data).  Furthermore, just putting symlinks there
# *also* seems complicated, because there's an apparent
# chicken-and-egg problem with RStudio trying to read data from the
# installed package location.... in order to install the symlinks.
#
# Another way is to check in the relatively small ReadsPerGene.out.tab
# files under 'data'.  However it would be good to see that the files
# have been correctly generated and yield expected output, first.  So,
# for now, the script must be manually updated to point to the local
# data location.

# dataDir <- "/Volumes/G-DRIVE\ mobile\ USB-C/Analysis/Testrun/STAR/results"
dataDir <- "/usr/local/DNFA-genfiles/data/r-extdata"

getFile <- function(stem) {
  name <- paste ('test', stem, 'ReadsPerGene.out.tab', sep = "")
  file.path(dataDir, name)
}

inputFiles <- lapply(stems, getFile)
test <- mergeFiles(inputFiles)
write.csv(test, file = 'data/gene-counts.csv')
