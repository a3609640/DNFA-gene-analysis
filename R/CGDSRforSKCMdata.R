# the following script generates the plot and statistics
# for DNFA expression and mutations data from recent TCGA pan-cancer dataset.
# install.packages("cgdsr")
library(car)
library(cgdsr)
library(cowplot)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(httr)
library(plyr)
library(reshape2)
library(stringr)
library(survival)
library(survMisc)
library(survminer)
library(tidyr)

  ### Create CGDS object
  ### mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
mycgds <- CGDS("http://www.cbioportal.org/")
  ### Get list of cancer studies at server
test(mycgds)
getCancerStudies(mycgds)

DNFA.gene <- c("ACLY", "ACSS2","ACACA", "SCD", "FASN", "ACSL1",
               "HMGCS1", "HMGCR", "MVK", "PMVK", "MITF", 
               "BRAF", "NRAS", "AKT1", "PTEN", "TP53")
names(DNFA.gene ) <- DNFA.gene

##############################################
## Get DNFA gene expression from SKCM group ##
##############################################
  ### skcm_case <- getCaseLists(mycgds, "skcm_tcga")
  ### skcm_tcga_all <- getCaseLists(mycgds, "skcm_tcga")[2, 1]
  ### Get available genetic profiles
  ### SKCMgeneticprofile <- getGeneticProfiles(mycgds, "skcm_tcga")
  ### getProfileData(mycgds,"genename","genetic profile IDs","A case list ID")
DNFA.RNAseq.data <- getProfileData(mycgds,
                                   DNFA.gene,
                                   "skcm_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
                                   "skcm_tcga_pan_can_atlas_2018_all")

################################################
## Get oncogene mutation data from SKCM group ##
################################################
  ### note there may be some internal bugs for the data labeled as NaN
  ### https://github.com/cBioPortal/cgdsr/issues/2
getmutations <- function(x) {
  mutations <- getProfileData(mycgds,
                              x,
                              "skcm_tcga_pan_can_atlas_2018_mutations",
                              "skcm_tcga_pan_can_atlas_2018_all")
  colnames(mutations) <- paste0(colnames(mutations), '.mutations')
  v <- rownames(mutations)
  ### each mutation column contains three types of data:
  ### mutation (V600E), NAN (wildtype), NA (not sequenced).
  relabel.mutations <- function(gene) {
    mutations[, gene] <- ifelse(
      mutations[, gene] == "NaN", "Wildtype", "Mutated")
  }
  ### use sapply , input as a matrix, and output as a matrix too.
  mutations <- sapply(colnames(mutations), relabel.mutations)
  mutations <- as.data.frame(mutations)
  ### the sapply function return a new matrix lacking row names.
  ### add row names back with the following two lines
  mutations2 <- cbind(Row.Names = v, mutations)
  ### mutations <- mutations2[ , -1]
  rownames(mutations) <- mutations2[ , 1]
  return(mutations)
}

  ### mutations.data <- getmutations("BRAF")
mutations.data <- getmutations(DNFA.gene)
  ### mutations.list <- c("BRAF", "NRAS", "AKT", "TP53")

###########################################
## Get oncogene CNV data from SKCM group ##
###########################################
getCNV <- function(x) {
  CNV <- getProfileData(mycgds,
                        x,
                        "skcm_tcga_pan_can_atlas_2018_gistic",
                        "skcm_tcga_pan_can_atlas_2018_all")
  v <- rownames(CNV)
  colnames(CNV) <- paste0(colnames(CNV), '.CNV')
  # change CNV column from numeric to factor
  factor.CNV <- function(gene) {
    CNV[, gene] <- as.factor(CNV[, gene])
  }
  # use sapply, input as a matrix, and output as a matrix too.
  CNV <- sapply(colnames(CNV), factor.CNV)
  CNV <- as.data.frame(CNV)
  CNV2 <- cbind(Row.Names = v, CNV)
  rownames(CNV) <- CNV2[ , 1]
  #  print(attributes(CNV))
  return(CNV)
}

  ### CNV.data <- getCNV("BRAF")
CNV.data <- getCNV(DNFA.gene)

###############################################################################
## examine the correlation between mutation status and DNFA expression level ##
###############################################################################
plot.mutations.RNAseq <- function(mutations, RNAseq) {
  mutations.DNFA.RNAseq <- cbind(mutations.data, DNFA.RNAseq.data)
  mutations.DNFA.RNAseq <- na.omit(mutations.DNFA.RNAseq)
  print(
    ggplot(mutations.DNFA.RNAseq,
           #use[ ,genemutations] not $genemutations for variable in a function
           aes(x     = mutations.DNFA.RNAseq[ ,mutations],
               y     = log2(mutations.DNFA.RNAseq[, RNAseq]),
               color = mutations.DNFA.RNAseq[ ,mutations])) +
      geom_boxplot(alpha = .01,
                   width = .5) +
      labs(x = mutations,
           y = paste("log2(", RNAseq, "RNA counts)")) +
      theme(axis.title        = element_text(face   = "bold",
                                             size   = 9,
                                             color  = "black"),
            axis.text         = element_text(size   = 9,
                                             hjust  = 1,
                                             face   = "bold",
                                             color  = "black"),
            axis.line.x       = element_line(color  = "black"),
            axis.line.y       = element_line(color  = "black"),
            panel.grid        = element_blank(),
            strip.text        = element_text(face   = "bold",
                                             size   = 9,
                                             colour = "black"),
            legend.position   = "none") +
      # stack the dots along the y-axis and group them along x-axis
      geom_dotplot(binaxis    = "y",
                   binwidth   = .1,
                   stackdir   = "center",
                   fill       = NA))
}

plot.mutations.RNAseq("NRAS.mutations", "SCD")
sapply(c("BRAF.mutations", "NRAS.mutations",
         "AKT1.mutations", "TP53.mutations"),
       function(x) mapply(plot.mutations.RNAseq,
                          x,
                          c("BRAF", "NRAS", "SCD")))
   ### or
sapply(c("BRAF", "NRAS", "SCD"),
       function(y) mapply(plot.mutations.RNAseq,
                          c("BRAF.mutations", "NRAS.mutations",
                            "AKT1.mutations", "TP53.mutations"),
                          y))


#######################################################################
##  check the data distribution, then choose the statistics methods  ##
#######################################################################
stats <- function(RNAseq, mutations) {
  mutations.DNFA.RNAseq <- cbind(mutations.data, DNFA.RNAseq.data)
  mutations.DNFA.RNAseq <- na.omit(mutations.DNFA.RNAseq)
  print(
    ggplot(mutations.DNFA.RNAseq,
           aes(x      = mutations.DNFA.RNAseq[, RNAseq],
               colour = mutations.DNFA.RNAseq[, mutations])) +
      geom_density() +
      scale_x_continuous(name = RNAseq) +
      scale_colour_discrete(name = mutations)) # change the legend title
  Mutated <- function(RNAseq, mutations) {
    mutations.DNFA.RNAseq[ ,RNAseq][
      mutations.DNFA.RNAseq[ ,mutations] == "Mutated"]
  }
  Wildtype <- function(RNAseq, mutations) {
    mutations.DNFA.RNAseq[ ,RNAseq][
      mutations.DNFA.RNAseq[ ,mutations] == "Wildtype"]
  }
  list <- list(Mutated(RNAseq, mutations),
               Wildtype(RNAseq, mutations))
  nameit <- function(RNAseq, mutations) {
    name <- c(paste(RNAseq, "with", mutations),
              paste(RNAseq, "without", mutations))}
  names(list) <- nameit(RNAseq, mutations)
  result <- lapply(list, shapiro.test)
  a <- t.test(
    mutations.DNFA.RNAseq[ ,RNAseq] ~ mutations.DNFA.RNAseq[ ,mutations]
  )
  # overwrite the data.name variable of t.test result
  a$data.name <- paste(RNAseq, 'expression by', mutations, 'status')
  b <- wilcox.test(
    mutations.DNFA.RNAseq[ ,RNAseq] ~ mutations.DNFA.RNAseq[ ,mutations]
  )
  b$data.name <- paste(RNAseq, 'expression by', mutations, 'status')
  print(result)
  print(a)
  print(b)
}

stats("SCD", "BRAF.mutations")
sapply(c("SCD", "FASN"),
       function(x)
         mapply(stats,
                x,
                c("BRAF.mutations",
                  "NRAS.mutations",
                  "AKT1.mutations",
                  "TP53.mutations")))
   ### or
sapply(c("BRAF.mutations",
         "NRAS.mutations",
         "AKT1.mutations",
         "TP53.mutations"),
       function(y)
         mapply(stats,
                c("SCD", "FASN"),
                y))

########################################################################
## compare the correlation between oncogene CNV and DNFA RNA-seq data ##
########################################################################
### make a large function to plot all genes
plot.CNV.RNAseq <- function(geneCNV, RNAseq) {
  CNV.DNFA.RNAseq <- cbind(CNV.data, DNFA.RNAseq.data)
  # na.omit cannot eleminate NaN here!
  toBeRemoved <- which(CNV.DNFA.RNAseq$BRAF.CNV == "NaN")
  CNV.DNFA.RNAseq <- CNV.DNFA.RNAseq[-toBeRemoved,]
  print(
    ggplot(CNV.DNFA.RNAseq,
           aes(x     = CNV.DNFA.RNAseq[, geneCNV],
               y     = log2(CNV.DNFA.RNAseq[, RNAseq]),
               color = CNV.DNFA.RNAseq[, geneCNV])) +
      geom_boxplot(alpha      = .01,
                   width      = .5) +
      geom_dotplot(binaxis    = "y",
                   binwidth   = .1,
                   stackdir   = "center",
                   fill       = NA) +
      scale_x_discrete(limits = c("-2", "-1", "0", "1", "2"), # skip NaN data
                       labels = c("Homdel", "Hetlos", "Diploid", "Gain", "Amp"),
                       drop   = FALSE) +
      labs(x = geneCNV,
           y = paste("log2(", RNAseq, "RNA counts)")) +
      theme(axis.title      = element_text(face   = "bold",
                                           size   = 9,
                                           color  = "black"),
            axis.text       = element_text(size   = 9,
                                           hjust  = 1,
                                           face   = "bold",
                                           color  = "black"),
            axis.line.x     = element_line(color  = "black"),
            axis.line.y     = element_line(color  = "black"),
            panel.grid      = element_blank(),
            strip.text      = element_text(face   = "bold",
                                           size   = 9,
                                           colour = "black"),
            legend.position = "none")
  )
  a <- leveneTest(CNV.DNFA.RNAseq[, RNAseq] ~ CNV.DNFA.RNAseq[, geneCNV])
  a$data.name <- paste(RNAseq,'expression by', geneCNV, 'status')
  b <- fligner.test(CNV.DNFA.RNAseq[, RNAseq] ~ CNV.DNFA.RNAseq[, geneCNV])
  b$data.name <- paste(RNAseq,'expression by', geneCNV, 'status')
  print(b)
  print(a)
}
### use all combinations of geneCNV and RNAseq
sapply(c("BRAF.CNV","NRAS.CNV", "PTEN.CNV"),
       function(x)
         mapply(plot.CNV.RNAseq, x, c("BRAF", "NRAS", "PTEN", "SCD", "FASN")))
### or
sapply(c("BRAF", "NRAS", "PTEN", "SCD", "FASN"),
       function(y)
         mapply(plot.CNV.RNAseq, c("BRAF.CNV","NRAS.CNV", "PTEN.CNV"), y))
### statistic comparision of SCD expression
### between PTEN heterdeletion and diploid
### (pten loss correlates with scd decrease?)

##################################################################
##  plot OS curve with clinic and mutation data from SKCM group ##
##################################################################
plot.mutation.SKCM.OS <- function(ge) {
  mycancerstudy <- getCancerStudies(mycgds)[
    grep("^skcm_tcga$", getCancerStudies(mycgds)$cancer_study_id), 1]
  ## ^ and $ indicate exact match
  mycaselist <- getCaseLists(mycgds, mycancerstudy)[4, 1]
  skcm.clinicaldata <- getClinicalData(mycgds, mycaselist)
  skcm.clinicaldata$rn <- rownames(skcm.clinicaldata)
  mutations.data <- na.omit(mutations.data)
  mutations.data$rn <- rownames(mutations.data)
  #  CNV.data <- getCNV(ge)
  # na.omit does not eleminate NaN in CNV.DNFA.RNAseq,
  # because NaN was treated as a factor in CNV.DNFA.RNAseq$BRAF.CNV
  #  CNV.data$rn <- rownames(CNV.data)
  #  CNVname <- paste0(ge, '.CNV')
  #  toBeRemoved <- which(CNV.data[[CNVname]] == "NaN")
  #  CNV.data <- CNV.data[-toBeRemoved,]
  df <- join_all(list(skcm.clinicaldata[c("OS_MONTHS",
                                          "OS_STATUS",
                                          "rn")],
                      mutations.data),
                 by   = "rn",
                 type = "full")
  df <- na.omit(df)
  df$SurvObj <- with(df, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))
  fit <- function(x) {
    survfit(SurvObj ~ df[[x]], data = df, conf.type = "log-log")
  }
  km <- fit(ge)
  black.bold.12pt <- element_text(face   = "bold",
                                  size   = 12,
                                  colour = "black")
  print(
    ggplot2::autoplot(km,
                      xlab = "Months",
                      ylab = "Survival Probability",
                      main = paste("Kaplan-Meier plot", ge)) +
      theme(axis.title      = black.bold.12pt,
            axis.text       = black.bold.12pt,
            axis.line.x     = element_line(color  = "black"),
            axis.line.y     = element_line(color  = "black"),
            panel.grid      = element_blank(),
            strip.text      = black.bold.12pt,
            legend.text     = black.bold.12pt ,
            legend.title    = black.bold.12pt ,
            legend.justification = c(1,1)) +
      scale_fill_discrete(name = ge))
  stats <- function(x) {
    survdiff(SurvObj ~ df[[x]], data = df, rho = 1)
  }
  print(ge)
  print(stats(ge))
}

plot.mutation.SKCM.OS("BRAF.mutations")
mutation.list <- c("BRAF.mutations",
                   "NRAS.mutations",
                   "AKT1.mutations",
                   "TP53.mutations")
names(mutation.list) <- mutation.list
sapply(mutation.list, plot.mutation.SKCM.OS)

#########################################################################
##  plot OS curve with clinic and DNFA expression data from SKCM group ##
#########################################################################
### PANCAN dataset from cBioportal does not offer clinical OS results  
### The following script uses OS data from TCGA provisional data
### and combine OS data with the expression data from pancan study
plot.DNFA.OS <- function(DNFA) {
  mycancerstudy <- getCancerStudies(mycgds)[
    grep("^skcm_tcga$", getCancerStudies(mycgds)$cancer_study_id), 1]
  mycaselist <- getCaseLists(mycgds, mycancerstudy)[4, 1] ### "skcm_tcga_all"
  skcm.clinicaldata <- getClinicalData(mycgds, mycaselist)
  skcm.clinicaldata$rn <- rownames(skcm.clinicaldata)
  skcm.RNAseq.data <- getProfileData(mycgds,
                                     DNFA,
                                     "skcm_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
                                     "skcm_tcga_pan_can_atlas_2018_all")
  skcm.RNAseq.data <- as.data.frame(skcm.RNAseq.data)
  skcm.RNAseq.data$rn <- rownames(skcm.RNAseq.data)
  df <- join_all(list(skcm.clinicaldata[c("OS_MONTHS", "OS_STATUS", "rn")],
                      skcm.RNAseq.data),
                 by   = "rn",
                 type = "full")
  df <- na.omit(df)
  df$Group[df[[DNFA]] < quantile(df[[DNFA]], prob = 0.2)] = "Bottom 20%"
  df$Group[df[[DNFA]] > quantile(df[[DNFA]], prob = 0.8)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))
  df <- na.omit(df)
  km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  black.bold.12pt <- element_text(face   = "bold",
                                  size   = 12,
                                  colour = "black")
  print(
    ggplot2::autoplot(km,
             xlab = "Months",
             ylab = "Survival Probability",
             main = paste("Kaplan-Meier plot", DNFA, "RNA expression")) +
      theme(axis.title           = black.bold.12pt,
            axis.text            = black.bold.12pt,
            axis.line.x          = element_line(color  = "black"),
            axis.line.y          = element_line(color  = "black"),
            panel.grid           = element_blank(),
            strip.text           = black.bold.12pt,
            legend.text          = black.bold.12pt ,
            legend.title         = black.bold.12pt ,
            legend.justification = c(1,1)))
  # rho = 1 the Gehan-Wilcoxon test
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 1)
  print(DNFA)
  print(stats)
}

plot.DNFA.OS("SCD")
sapply(DNFA.gene, plot.DNFA.OS)



