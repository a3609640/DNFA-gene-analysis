
# Draw boxplots to compare the expression of interested genes from 
# different tumor types in the TCGA database.

# Combine the individual gene expression together into a big file
# to compare the correlation between indicated genes across
# different tumor types in TCGA 

```{r, echo=TRUE}
library(stringr)
library(reshape2)
library(ggplot2)
```

```{r, echo=TRUE}
dataDir <- "~/TCGA cBioportal"
setwd("~/TCGA cBioportal")
```

# Extract the name of a gene from that gene's data file.
```{r}
getGeneName <- function(geneDataFile) {
  # The gene name is part of the data file name, so we just
  # discard the parts of the file name that we don't need.
  geneDataFile <- str_replace_all(geneDataFile, dataDir, "");
  geneDataFile <- str_replace_all(geneDataFile, ".csv", "");
  geneDataFile <- str_replace_all(geneDataFile, "/", "");
  return(geneDataFile);
}
```

```{r}

```
# This function does the following:
# 1. Given the name of a CSV file for a particular gene,
#    read a data table from that file.
# 2. Modify the table by removing columns 3 and 4.
# 3. Modifies the table's "cancer study" column to trim some
#    redundant phrases.
# 4. Return the modified data table (with removed
#    and renamed columns) to the caller.

```{r}
getGeneTableFromFile <- function(geneDataFile) {
  geneDataTable <- read.csv(geneDataFile)
  geneDataTable <- geneDataTable[, -3:-4]
  #  remove (TCGA, Provisional) from the string, \\ to remove special character such as (), 
  geneDataTable$Cancer.Study <- str_replace_all(geneDataTable$Cancer.Study, "\\(TCGA, Provisional\\)", "")
  return(geneDataTable)
}
```

# This function does the following:
# 1. Given the name of a CSV file containing data for a gene,
#    gets a corresponding data table as provided by changeName().
# 2. Re-orders the rows of the table, sorting the cancer groups
#    by expression of the gene, from lowest to highest.
# 3. Generate a box plot of the re-ordered data, supplying
#    axis and legend labels as appropriate.

```{r}
plotGene <- function(geneDataFile) {
  geneName <- getGeneName(geneDataFile)
  geneDataTable <- getGeneTableFromFile(geneDataFile)
  # reorder the data according to the mean of RNA-seq value
  mean <- within(geneDataTable, Cancer.Study <-  reorder(Cancer.Study, log2(Value), median))
  # fill=Cancer.Study color each boxplot with different color
  ggplot(mean, aes(x=Cancer.Study, y=log2(Value), fill=Cancer.Study)) + 
    geom_violin(alpha = .5, trim=FALSE)+
    geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))+
    # scale_fill_manual(values = great_color) +
    # geom_dotplot(binaxis="y", binwidth=.05,stackdir="centerwhole")+ 
    theme_bw()+
    labs(x = "Tumor types (TCGA)",
         y = paste("log2(", geneName, "RNA counts)")) +
    theme(axis.title=element_text(face="bold",size=9,color="black"),
          axis.text=element_text(size=9,angle = 45, hjust = 1, face="bold",color="black"),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 9, colour = "black"),
          legend.position = "none")}
```

plotGene(geneDataFile="MITF.csv")

```{r, include = TRUE, echo = TRUE}
gene.list <- list.files(dataDir, '.csv')
lapply(gene.list, plotGene)
```

# This function creates a dataset, and merges data
# from each file in fileList into it.  The data are
# joined by row names.  The 'Value' column in each
# merged data set is renamed with the name of its
# gene.

```{r}
mergeGeneDataFiles <- function(fileList) {

  # TODO(su) find the right way to initialize 'dataset'
  dataset <- NULL

  # TODO(su) Current bug:
  # --
  # Error in data.frame(..., check.names = FALSE) : 
  #   arguments imply differing number of rows: 0, 9721
  # --
  # Need to find a way to merge data sets with differing
  # numbers of rows, OR modify getGeneTableFromFile so that
  # it returns tables with uniform sizes.
  for (geneDataFile in gene.list){
    geneName <- getGeneName(geneDataFile)
    geneData <- getGeneTableFromFile(geneDataFile)
    colnames(geneData) <- str_replace_all(colnames(geneData), "Value", geneName)
    # library(plyr)
    # dataset <- rbind.fill(dataset, geneData)
    # dataset<-cbind(dataset, geneData)
  }
  return(dataset)
}
```

if(interactive()) {
## Use list.files to obtain all CSV files in 'dataDir'
## (these were downloaded from cbioportal website),
## then generate box plots for all of them.
gene.list <- list.files(dataDir, '.csv')
lapply(gene.list, plotGene)

mergedDataSet <- mergeGeneDataFiles(gene.list)

### Below this line are some code snippets being used for manual
### inspection of the data.

cor.test(mergedDataSet$FASN,mergedDataSet$SCD, method = "pearson")
cor.test(mergedDataSet$FASN,mergedDataSet$SREBF1, method = "pearson")
cor.test(mergedDataSet$MITF,mergedDataSet$SCD, method = "pearson")
p <- ggplot(mergedDataSet, aes(FASN, SCD))  
p + geom_point()
p + geom_point(aes(colour= factor(Cancer.Study)))+ facet_wrap(~Cancer.Study)  

dataset.skin <- subset(dataset, Cancer.Study=='Skin Cutaneous Melanoma (TCGA, Provisional)')
cor.test(mergedDataSet.skin$FASN,mergedDataSet.skin$SCD, method = "pearson")
cor.test(mergedDataSet.skin$FASN,mergedDataSet.skin$SREBF1, method = "pearson")
cor.test(mergedDataSet.skin$MITF,mergedDataSet.skin$SCD, method = "pearson")

p <- ggplot(mergedDataSet.skin, aes(log2(FASN), log2(SCD)))  
p + geom_point()

p <- ggplot(mergedDataSet.skin, aes(log2(MITF), log2(SCD)))  
p + geom_point()


}