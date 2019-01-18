# the following script generates the plot and statistics
# for DNFA expression and mutations data from TCGA provisional dataset.
# install.packages("cgdsr")
library(car)
library(cgdsr)
library(ggfortify)
library(ggplot2)
library(httr)
library(plyr)
library(reshape2)
library(stringr)
library(survival)

# Create CGDS object
mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
# mycgds = CGDS("http://www.cbioportal.org/")
test(mycgds)
# Get list of cancer studies at server
getCancerStudies(mycgds)
# Get cases from TCGA provisional studies only
DNFA.gene <- c("ACACA", "SCD", "ACLY", "FASN", "ACSS2", "MITF")
names(DNFA.gene) <- DNFA.gene

#############################################################
## plot DNFA RNASeq data from all TCGA cancer study groups ##
#############################################################
plot.DNFA.tcga <- function(DNFA){
  # Get DNFA RNAseq data from all TCGA study groups
  tcga.pro.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)", getCancerStudies(mycgds)$name), ]
  # "tcag_study_list" contains all the tcga cancer studies
  tcga.study.list <- tcga.pro.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  # use lappy to pull out all the caselists within tcga.study.list
  # because we named each elements in tcga.study.list,
  # lappy will return a large list, each element (with a cancer study name)
  # in that list is a data-table
  tcga.pro.caselist <- lapply(tcga.study.list, caselist)
  tcga.pro.geneticprofile <- lapply(tcga.study.list, geneticprofile)
  # for example, tcga.pro.caselist[[1]] shows the dataframe of caselist
  # in laml study group.
  # to choose case_list_id that is labeled with laml_tcga_rna_seq_v2_mrna,
  # we use the following tcag_provisional_caselist[[1][8,1]
  # a <- tcga.pro.caselist[[1]][
  # grep("tcga_rna_seq_v2_mrna", tcga.pro.caselist[[1]]$case_list_id),
  # ][1,1]
  # b <- tcga.pro.geneticprofile[[1]][
  # grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
  # tcga.pro.geneticprofile[[1]]$genetic_profile_name), ][1,1]
  # how do we do this for all study groups from [[1]] to  [[32]]?
  caselist.RNAseq <- function(x) {
    tcga.pro.caselist[[x]][
      grep("tcga_rna_seq_v2_mrna",
           tcga.pro.caselist[[x]]$case_list_id), ][1, 1]
    }
  geneticprofile.RNAseq <- function(x) {
    tcga.pro.geneticprofile[[x]][
  # double backslash \\ suppress the special meaning of ( )
  # in regular expression
      grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
           tcga.pro.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
    }
  # test the functions: caselist.RNAseq () and geneticprofile.RNAseq ()
  # caselist.RNAseq = caselist.RNAseq ('acc_tcga')
  # geneticprofile.RNAseq = geneticprofile.RNAseq ('acc_tcga')
  # Wrap two functions: geneticprofile.RNAseq(x), caselist.RNAseq(x)
  # within TCGA_ProfileData_RNAseq(x)
  tcga.profiledata.RNAseq <- function(genename, geneticprofile, caselist) {
    getProfileData(mycgds,
                   genename,
                   geneticprofile,
                   caselist)
    }
  DNFA.tcga.RNAseq <- function(x, y) {
    tcga.profiledata.RNAseq(x,
                            geneticprofile.RNAseq(y),
                            caselist.RNAseq(y))
    }
  DNFA.RNAseq.tcga.all <- function(x) {
    test <- lapply(tcga.study.list,
                   function(y) mapply(DNFA.tcga.RNAseq, x, y))
    df2 <- melt(test)
    colnames(df2) <- c("RNAseq", "DNFAgene", "TCGAstudy")
    df2 <- data.frame(df2)
    }
  df2 <- DNFA.RNAseq.tcga.all(DNFA)
  df2$DNFAgene <- as.factor(df2$DNFAgene)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  # plot DNFA gene expression across all TCGA groups ##
  m <- paste0(DNFA, ".", DNFA)
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
           y = paste0("log2(", DNFA, " RNA counts)")) +
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

plot.DNFA.tcga("SCD")
sapply(DNFA.gene, plot.DNFA.tcga)

##############################################
## Get DNFA gene expression from SKCM group ##
##############################################
# skcm_case <- getCaseLists(mycgds, "skcm_tcga")
# skcm_tcga_all <- getCaseLists(mycgds, "skcm_tcga")[2, 1]
# Get available genetic profiles
# SKCMgeneticprofile <- getGeneticProfiles(mycgds, "skcm_tcga")
DNFA.RNAseq.data <- getProfileData(mycgds,
                                   c("ACACA", "FASN", "SCD", "ACLY", "ACSS2",
                                     "HMGCS1", "HMGCR", "SREBF1", "SREBF2",
                                     "MITF", "BRAF", "NRAS", "PTEN"),
                                   "skcm_tcga_rna_seq_v2_mrna",
                                   "skcm_tcga_all")

################################################
## Get oncogene mutation data from SKCM group ##
################################################
# note there may be some internal bugs for the data labeled as NaN
# https://github.com/cBioPortal/cgdsr/issues/2
getmutations <- function() {
  mutations <- getProfileData(mycgds,
                              c("BRAF", "NRAS", "AKT", "TP53"),
                              "skcm_tcga_mutations",
                              "skcm_tcga_all")
  colnames(mutations) <- paste0(colnames(mutations), '.mutations')
  v <- rownames(mutations)
  # each mutation column contains three types of data:
  # mutation (V600E), NAN (wildtype), NA (not sequenced).
  relabel.mutations <- function(gene) {
    mutations[, gene] <- ifelse(
      mutations[, gene] == "NaN", "Wildtype", "Mutated")
    }
  # use sapply , input as a matrix, and output as a matrix too.
  mutations <- sapply(colnames(mutations), relabel.mutations)
  mutations <- as.data.frame(mutations)
  # the sapply function return a new matrix lacking row names.
  # add row names back with the following two lines
  mutations2 <- cbind(Row.Names = v, mutations)
  # mutations <- mutations2[ , -1]
  rownames(mutations) <- mutations2[ , 1]
  return(mutations)
  }

# mutations.data <- getmutations("BRAF")
mutations.data <- getmutations()
# mutations.list <- c("BRAF", "NRAS", "AKT", "TP53")

###########################################
## Get oncogene CNV data from SKCM group ##
###########################################
getCNV <- function(x) {
  CNV <- getProfileData(mycgds,
                        x,
                        "skcm_tcga_gistic",
                        "skcm_tcga_all")
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

# CNV.data <- getCNV("BRAF")
CNV.data <- getCNV(c("BRAF", "NRAS", "PTEN"))

###########################################################################
## examine the correlation between mutation status and DNFA RNASeq level ##
###########################################################################
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
sapply(c("BRAF.mutations","NRAS.mutations",
         "AKT1.mutations", "TP53.mutations"),
       function(x) mapply(plot.mutations.RNAseq,
                          x,
                          c("BRAF", "NRAS", "SCD")))
# or
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
# or
sapply(c("BRAF.mutations",
         "NRAS.mutations",
         "AKT1.mutations",
         "TP53.mutations"),
       function(y)
         mapply(stats,
                c("SCD", "FASN"),
                y))

#############################################################################
## compare the correlation between oncogene CNV and DNFA gene RNA-seq data ##
#############################################################################
# make a large function to plot all genes
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
# use all combinations of geneCNV and RNAseq
sapply(c("BRAF.CNV","NRAS.CNV", "PTEN.CNV"),
       function(x)
         mapply(plot.CNV.RNAseq, x, c("BRAF", "NRAS", "PTEN", "SCD", "FASN")))
# or
sapply(c("BRAF", "NRAS", "PTEN", "SCD", "FASN"),
       function(y)
         mapply(plot.CNV.RNAseq, c("BRAF.CNV","NRAS.CNV", "PTEN.CNV"), y))
# statistic comparision of SCD expression
# between PTEN heterdeletion and diploid
# (pten loss correlates with scd decrease?)

##################################################################
##  Kaplan-Meier curve with clinic and mutation data from SKCM  ##
##################################################################
plot.km.mut.skcm <- function(ge) {
  mycancerstudy <- getCancerStudies(mycgds)[194, 1]
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
    autoplot(km,
             xlab = "Months",
             ylab = "Survival Probability",
             main = "Kaplan-Meier plot") +
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

plot.km.mut.skcm("BRAF.mutations")
mutation.list <- c("BRAF.mutations",
                   "NRAS.mutations",
                   "AKT1.mutations",
                   "TP53.mutations")
names(mutation.list) <- mutation.list
sapply(mutation.list, plot.km.mut.skcm)

#####################################################################
##  Kaplan-Meier curve with clinic and DNFA RNASeq data from SKCM  ##
#####################################################################
plot.km.DNFA.skcm <- function(DNFA) {
  mycancerstudy <- getCancerStudies(mycgds)[196, 1]        # "skcm_tcga"
  mycaselist <- getCaseLists(mycgds, mycancerstudy)[4, 1]  # "skcm_tcga_all"
  skcm.clinicaldata <- getClinicalData(mycgds, mycaselist)
  skcm.clinicaldata$rn <- rownames(skcm.clinicaldata)
  skcm.RNAseq.data <- getProfileData(mycgds,
                                     DNFA,
                                     "skcm_tcga_rna_seq_v2_mrna",
                                     "skcm_tcga_all")
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
    autoplot(km,
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

plot.km.DNFA.skcm("ACACA")
sapply(DNFA.gene, plot.km.DNFA.skcm)


##############################################################################
## Kaplan-Meier curve with clinic and DNFA RNAseq data from all TCGA groups ##
##############################################################################
plot.km.all.tcga <- function(DNFA) {
#  mycgds <- CGDS("http://www.cbioportal.org/")
#  test(mycgds)
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  tcga.pro.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)",
         getCancerStudies(mycgds)$name), ]
  tcga.study.list <- tcga.pro.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  tcga.pro.caselist <- lapply(tcga.study.list,
                              caselist)
  tcga.pro.geneticprofile <- lapply(tcga.study.list,
                                    geneticprofile)
  caselist.RNAseq <- function(x) {
    tcga.pro.caselist[[x]][grep("tcga_rna_seq_v2_mrna",
                                tcga.pro.caselist[[x]]$case_list_id), ][1, 1]
    }
  geneticprofile.RNAseq <- function(x) {
    tcga.pro.geneticprofile[[x]][
      grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
           tcga.pro.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
    }
  tcga.profiledata.RNAseq <- function(genename,
                                      geneticprofile,
                                      caselist) {
    getProfileData(mycgds,
                   genename,
                   geneticprofile,
                   caselist)
    }
  tcga.gene.RNAseq <- function(x, y) {
    tcga.profiledata.RNAseq(x,
                            geneticprofile.RNAseq(y),
                            caselist.RNAseq(y))
  }
  tcga.DNFA.RNAseq <- function(y) {
    tcga.gene.RNAseq(x = DNFA, y)
    }
## test1 <- SCD.tcga.RNAseq(y = "skcm_tcga")
## test ## try to keep patient ID in the rowname
  test2 <- lapply(tcga.study.list, tcga.DNFA.RNAseq)
  for(x in 1:32)
    {
    test2[[x]]$case.id <- rownames(test2[[x]])
    message("test2 = ", x)
    }
  df1 <- melt(test2)
  colnames(df1) <- c("case.id",
                     "DNFAgene",
                     "RNAseq",
                     "TCGAstudy")
  all.tcga.DNFA.RNAseq <- data.frame(df1)
  all.tcga.DNFA.RNAseq$TCGAstudy <- as.factor(all.tcga.DNFA.RNAseq$TCGAstudy)
  message("RNAseq data retrieved")
  ##### retrieve clinic data from all tcga groups #####
  tcga.clinic.data <- function(x) {
    print(x)
    url <- function(x){
      url <- "http://www.cbioportal.org/webservice.do?cmd=getClinicalData&case_set_id="
      url <- paste0(url, x, "_all")
      return(url)
      }
# testurl <- url("acc_tcga")
# tesereq <- GET(url("acc_tcga"))
  req <- function(x) {GET(url(x))}
# req <- req("acc_tcga")
  clinical_data <- function(x) {content(req(x),
                                        type      = 'text/tab-separated-values',
                                        col_names = T,
                                        col_types = NULL)}
  data <- clinical_data(x)
  data <- data[c("OS_MONTHS",
                 "OS_STATUS",
                 "CASE_ID")]
  }
# three datasets donot have OS data and cause bugs remove them
  bug.data.set <- names(tcga.study.list) %in% c("meso_tcga", "pcpg_tcga", "ucs_tcga")
  tcga.study.list <- tcga.study.list[!bug.data.set]
  all.tcga.clinic.data <- lapply(tcga.study.list, tcga.clinic.data)
  all.tcga.clinic.data <- melt(all.tcga.clinic.data)
  all.tcga.clinic.data <- all.tcga.clinic.data[c("OS_STATUS",
                                                 "CASE_ID",
                                                 "value",
                                                 "L1")]
  colnames(all.tcga.clinic.data) <- c("OS_STATUS",
                                      "case.id",
                                      "OS_MONTHS",
                                      "TCGAstudy")
  all.tcga.clinic.data$case.id <- str_replace_all(all.tcga.clinic.data$case.id,
                                                  '-',
                                                  '.')
  message("clinical data retrieved")
  df <- join_all(list(all.tcga.clinic.data[c("OS_MONTHS",
                                             "OS_STATUS",
                                             "case.id")],
                      all.tcga.DNFA.RNAseq[c("case.id",
                                             "RNAseq",
                                             "TCGAstudy")]),
                 by   = "case.id",
                 type = "full")
  df <- na.omit(df)
  message("clinical and RNAseq data combined")
  df$Group[df$RNAseq < quantile(df$RNAseq, prob = 0.2)] = "Bottom 20%"
  df$Group[df$RNAseq > quantile(df$RNAseq, prob = 0.8)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))
#  df <- na.omit(df)
  km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  black.bold.12pt <- element_text(face   = "bold",
                                  size   = 12,
                                  colour = "black")
  print(
    autoplot(km,
             xlab = "Months",
             ylab = "Survival Probability",
             main = paste("Kaplan-Meier plot", DNFA, "RNA expression"),
             xlim = c(0, 250)) +
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

plot.km.all.tcga("SCD")
sapply(DNFA.gene, plot.km.all.tcga)
