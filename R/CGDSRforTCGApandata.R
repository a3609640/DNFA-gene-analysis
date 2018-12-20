# the following script generates the plot and statistics
# for DNFA expression and mutations data from TCGA dataset.
# install.packages("cgdsr")
library(car)
library(cgdsr)
library(ggfortify)
library(ggplot2)
library(plyr)
library(reshape2)
library(survival)
library(UCSCXenaTools)

# Create CGDS object
mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
test(mycgds)
# Get list of cancer studies at server
getCancerStudies(mycgds)
# Get cases from TCGA provisional studies only
#####################################################################
## Get DNFA gene expression data from all TCGA cancer study groups ##
#####################################################################
getDNFAdata <- function(ge) {
  tcga.pan.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
  # "tcag_study_list" is a vector containing all the tcga cancer studies
  # that I would to analyze for DNFA gene expression
  tcga.study.list <- tcga.pan.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  # use lappy to pull out all the caselists within tcga.study.list
  # because we named each elements in tcga.study.list
  # (names(tcga.study.list) <- tcga.study.list),
  # lappy will return a large list, each element (with a cancer study name)
  # in that list is a data-table
  tcga.pan.caselist <- lapply(tcga.study.list, caselist)
  tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  # for example, tcga.pro.caselist[[1]] shows the dataframe of caselist
  # in laml study group.
  # to choose case_list_id that is labeled with laml_tcga_rna_seq_v2_mrna,
  # we use the following tcag_provisional_caselist[[1][8,1]
  # a <- tcga.pro.caselist[[1]][
  # grep("tcga_rna_seq_v2_mrna", tcag_provisional_caselist[[1]]$case_list_id),
  # ][1,1]
  # b <- tcga.pro.geneticprofile[[1]][
  # grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
  # tcga.pro.geneticprofile[[1]]$genetic_profile_name), ][1,1]
  # how do we do this for all study groups from [[1]] to  [[32]]?
  caselist.RNAseq <- function(x) {
    tcga.pan.caselist[[x]][
      grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
           tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.RNAseq <- function(x) {
    tcga.pan.geneticprofile[[x]][
      # double backslash \\ suppress the special meaning of ( )
      # in regular expression
      grep("mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)",
           tcga.pan.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
  # test the functions: caselist.RNAseq () and geneticprofile.RNAseq ()
  # caselist.RNAseq = caselist.RNAseq ('acc_tcga')
  # geneticprofile.RNAseq = geneticprofile.RNAseq ('acc_tcga')
  # We wrap two functions: geneticprofile.RNAseq(x), caselist.RNAseq(x)
  # within TCGA_ProfileData_RNAseq(x)
  tcga.profiledata.RNAseq <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds, genename, geneticprofile, caselist)
  }
  DNFA.tcga.RNAseq <- function(x, y) {
    tcga.profiledata.RNAseq(x, geneticprofile.RNAseq(y), caselist.RNAseq(y))
  }
  DNFA.RNAseq.all.tcga.studies <- function(x) {
    test <-
      lapply(tcga.study.list, function(y) mapply(DNFA.tcga.RNAseq, x, y))
    df2 <- melt(test)
    colnames(df2) <- c("RNAseq", "DNFAgene", "TCGAstudy")
    as.factor(df2$L2)
    as.factor(df2$L1)
    df2 <- data.frame(df2)
  }
  df2 <- DNFA.RNAseq.all.tcga.studies(ge)
  return(df2)
}

DNFA.gene <- c("SCD", "FASN", "ACLY", "ACSS2")
names(DNFA.gene ) <- DNFA.gene
df2 <- getDNFAdata(DNFA.gene)

######################################################
## plot DNFA gene expression across all TCGA groups ##
######################################################
plotDNFA <- function(x) {
  m <- paste0(x, ".", x)
  mean <- within(df2[df2$DNFAgene == m,], # TCGAstudy is one column in df2
                 TCGAstudy <- reorder(TCGAstudy, log2(RNAseq), median))
  print(
    ggplot(mean,
           aes(x     = TCGAstudy,
               y     = log2(RNAseq),
               color = TCGAstudy)) +
      geom_boxplot(alpha    = .01,
                   width    = .5,
                   position = position_dodge(width = .9)) +
      labs(x = "Tumor types (TCGA)",
           y = paste0("log2(", x, " RNA counts)")) +
      theme(axis.title  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            axis.text.x = element_text(size   = 9,
                                       angle  = 45,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.text.y = element_text(size   = 9,
                                       angle  = 0,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.line.x = element_line(color  = "black"),
            axis.line.y = element_line(color  = "black"),
            panel.grid  = element_blank(),
            strip.text  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            legend.position = "none"))
}

plotDNFA("SCD")
sapply(DNFA.gene, plotDNFA)