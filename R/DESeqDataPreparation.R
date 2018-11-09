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
processFile <- function(fileBaseDir, fileName) {
  fullFileName <- paste(fileBaseDir, fileName, sep='/')
  table <- read.table(fullFileName, stringsAsFactors=T)
  
  # Reframe the data, taking row names from first column.
  table.with.rownames <- data.frame(table[,-1], row.names=table[,1])
  
  # Generate 'test6.4.' from 'test6-4ReadsPerGene.out.tab'.
  newColBaseName <- paste(substring(fileName, 1, 5),
                          substring(fileName, 7, 7),
                          '',  # gives us a trailing '.'
                          sep='.')
  
  # Rename the columns as test6.4.1, test6.4.2, test6.4.3.
  colnames(table.with.rownames) <- paste0(newColBaseName, 1:3)
  return(table.with.rownames)
}

# Merge reframed data from all input files, and add 'symbol' column.
mergeFiles <- function(baseDir, files) {
  processedData <- lapply(inputFiles, processFile, fileBaseDir=baseDir)
  test <- do.call(cbind.data.frame, processedData)
  
  test$symbol <- mapIds(org.Hs.eg.db,
                        keys=row.names(test),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
  return(test)
}

dataDir <- "/Volumes/G-DRIVE\ mobile\ USB-C/Analysis/Testrun/STAR/results"

inputFiles<-c(
  "test1-1ReadsPerGene.out.tab",
  "test1-2ReadsPerGene.out.tab",
  "test1-3ReadsPerGene.out.tab",
  "test1-4ReadsPerGene.out.tab",
  "test2-1ReadsPerGene.out.tab",
  "test2-2ReadsPerGene.out.tab",
  "test2-3ReadsPerGene.out.tab",
  "test2-4ReadsPerGene.out.tab",
  "test3-1ReadsPerGene.out.tab",
  "test3-2ReadsPerGene.out.tab",
  "test3-3ReadsPerGene.out.tab",
  "test3-4ReadsPerGene.out.tab",
  "test4-1ReadsPerGene.out.tab",
  "test4-2ReadsPerGene.out.tab",
  "test4-3ReadsPerGene.out.tab",
  "test4-4ReadsPerGene.out.tab",
  "test5-1ReadsPerGene.out.tab",
  "test5-2ReadsPerGene.out.tab",
  "test5-3ReadsPerGene.out.tab",
  "test5-4ReadsPerGene.out.tab",
  "test6-1ReadsPerGene.out.tab",
  "test6-2ReadsPerGene.out.tab",
  "test6-3ReadsPerGene.out.tab",
  "test6-4ReadsPerGene.out.tab"
)

test <- mergeFiles(dataDir, inputFiles)
