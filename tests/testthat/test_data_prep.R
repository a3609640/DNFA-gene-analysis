library(base)
library(testthat)
library(DNFAGeneAnalysis)

setup({
  #  currently just a placeholder
})
teardown({
  #  currently just a placeholder
})

test_that("base column name generation works", {
  expect_equal("test1.1.", makeColumnBaseNameFromFileName("test1_S2_L001ReadsPerGene.out.tab"))
  expect_equal("test1.2.", makeColumnBaseNameFromFileName("test1_S2_L002ReadsPerGene.out.tab"))

  # TODO(dlroxe) Implementation isn't yet fancy enough for double-digit test numbers.
  # expect_equal("test1.12.", makeColumnBaseNameFromFileName("test1-12ReadsPerGene.out.tab"))
})

test_that("processFile works", {
   file1 <- normalizePath(file.path("..", "testdata", "test1_S2_L001ReadsPerGene.out.tab"))
   expect_true(file.exists(file1))
   processedData <- processReadsFile(file1)

   dimensions <- dim(processedData)
   rowCount <- nrow(processedData)
   size <- length(processedData)
   colNames <- colnames(processedData)
   
   expect_equal(c(7,1), dimensions)
   expect_equal(7, rowCount)
   expect_equal(1, size)
   expect_equal(c('test1.1.3'), colNames)
})

test_that("mergeFiles works", {
   file1 <- file.path("..", "testdata", "test1_S2_L001ReadsPerGene.out.tab")
   file2 <- file.path("..", "testdata", "test1_S2_L002ReadsPerGene.out.tab")
   expect_true(file.exists(file1))
   expect_true(file.exists(file2))
   mergedData <- mergeFiles(c(file1, file2))

   print(mergedData)
   dimensions <- dim(mergedData)
   rowCount <- nrow(mergedData)
   size <- length(mergedData)
   colNames <- colnames(mergedData)
 
   expect_equal(c(7,3), dimensions)
   expect_equal(7, rowCount)
   expect_equal(3, size)
   expect_equal(c('test1.1.3', 'test1.2.3', 'symbol'), colNames)
})

