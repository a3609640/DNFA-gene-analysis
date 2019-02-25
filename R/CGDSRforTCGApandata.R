# the following script generates the plot and statistics
# for DNFA expression and mutations data from recent TCGA pan-cancer dataset.
# install.packages("cgdsr")
library(car)
library(cgdsr)
library(ggfortify)
library(ggplot2)
library(plyr)
library(reshape2)
library(survival)
library(survminer)

# Create CGDS object
mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
mycgds = CGDS("http://www.cbioportal.org/")
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
  # "tcag_study_list" contains all the tcga cancer studies
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
  # for example, tcga.pan.caselist[[1]] shows the dataframe of caselist
  # in laml study group.
  # to choose case_list_id that is labeled with laml_tcga_rna_seq_v2_mrna,
  # we use the following tcag_provisional_caselist[[1][8,1]
  # a <- tcga.pan.caselist[[1]][
  # grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
  #      tcga.pan.caselist[[1]]$case_list_id), ][1,1]

  # b <- tcga.pan.geneticprofile[[1]][      # laml_tcga_pan_can_atlas_2018
  # grep("mRNA Expression, RSEM",
  #      tcga.pan.geneticprofile[[1]]$genetic_profile_name), ][1,1]
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
      grep("mRNA Expression, RSEM",
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
    test <- lapply(tcga.study.list, function(y) mapply(DNFA.tcga.RNAseq, x, y))
    df2 <- melt(test)
    colnames(df2) <- c("RNAseq", "DNFAgene", "TCGAstudy")
    as.factor(df2$L2)
    as.factor(df2$L1)
    df2 <- data.frame(df2)
    }
  df2 <- DNFA.RNAseq.all.tcga.studies(ge)
  return(df2)
}

DNFA.gene <- c("ACLY", "ACSS2","ACACA", "SCD", "FASN", "ACSL1", 
               "HMGCS1", "HMGCR", "MVK")
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

#############################################################
## plot DNFA RNASeq data from pan TCGA cancer study groups ##
#############################################################
plot.DNFA.pan.tcga <- function(EIF){
  # Get EIF RNAseq data from all TCGA study groups
  tcga.pan.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)", getCancerStudies(mycgds)$name), ]
  # "tcag_study_list" contains all the tcga cancer studies
  tcga.study.list <- tcga.pan.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  # use lappy to pull out all the caselists within tcga.study.list
  # because we named each elements in tcga.study.list,
  # lappy will return a large list, each element (with a cancer study name)
  # in that list is a data-table
  tcga.pan.caselist <- lapply(tcga.study.list, caselist)
  tcga.pan.geneticprofile <- lapply(tcga.study.list, geneticprofile)
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
    tcga.pan.caselist[[x]][
      grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
           tcga.pan.caselist[[x]]$case_list_id), ][1, 1]
  }
  geneticprofile.RNAseq <- function(x) {
    tcga.pan.geneticprofile[[x]][
      # double backslash \\ suppress the special meaning of ( )
      # in regular expression
      grep("mRNA Expression, RSEM",
           tcga.pan.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
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
  EIF.tcga.RNAseq <- function(x, y) {
    tcga.profiledata.RNAseq(x,
                            geneticprofile.RNAseq(y),
                            caselist.RNAseq(y))
  }
  EIF.RNAseq.tcga.all <- function(x) {
    test <- lapply(tcga.study.list,
                   function(y) mapply(EIF.tcga.RNAseq, x, y))
    df2 <- melt(test)
    colnames(df2) <- c("RNAseq", "EIFgene", "TCGAstudy")
    df2 <- data.frame(df2)
  }
  df2 <- EIF.RNAseq.tcga.all(EIF)
  df2$EIFgene <- as.factor(df2$EIFgene)
  df2$TCGAstudy <- as.factor(df2$TCGAstudy)
  df2 <- na.omit(df2)
  # plot EIF gene expression across all TCGA groups ##
  m <- paste0(EIF, ".", EIF)
  mean <- within(df2[df2$EIFgene == m,], # TCGAstudy is one column in df2
                 TCGAstudy <- reorder(TCGAstudy, log2(RNAseq), median))
  a <- levels(mean$TCGAstudy)
  # with highlight on skin cancer
  colors <- ifelse(a == "skcm_tcga_pan_can_atlas_2018", "red", "black")
  print(
    ggplot(mean,
           aes(x        = TCGAstudy,
               y        = log2(RNAseq),
               color    = TCGAstudy)) +
      geom_boxplot(alpha    = .01,
                   width    = .5,
                   position = position_dodge(width = .9)) +
      coord_flip() +
      labs(x = "Tumor types (TCGA)",
           y = paste0("log2(", EIF, " RNA counts)")) +
      theme(axis.title  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            axis.text.x = element_text(size   = 9,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = "black"),
            axis.text.y = element_text(size   = 9,
                                       angle  = 0,
                                       hjust  = 1, # 1 means right-justified
                                       face   = "bold",
                                       color  = colors),
            axis.line.x = element_line(color  = "black"),
            axis.line.y = element_line(color  = "black"),
            panel.grid  = element_blank(),
            strip.text  = element_text(face   = "bold",
                                       size   = 9,
                                       color  = "black"),
            legend.position = "none"))
}

plot.DNFA.pan.tcga("HMGCR")
sapply(DNFA.gene, plot.DNFA.pan.tcga)

##############################################
## Get DNFA gene expression from SKCM group ##
##############################################
# skcm_case <- getCaseLists(mycgds, "skcm_tcga")
# skcm_tcga_all <- getCaseLists(mycgds, "skcm_tcga")[2, 1]
# Get available genetic profiles
# SKCMgeneticprofile <- getGeneticProfiles(mycgds, "skcm_tcga")
# getProfileData(mycgds,"genename","genetic profile IDs","A case list ID")
DNFA.RNAseq.data <- getProfileData(mycgds,
                                   c("ACACA", "FASN", "SCD", "ACLY", "ACSS2",
                                     "HMGCS1", "HMGCR", "SREBF1", "SREBF2",
                                     "MITF", "BRAF", "NRAS", "PTEN"),
                                   "skcm_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
                                   "skcm_tcga_pan_can_atlas_2018_all")

################################################
## Get oncogene mutation data from SKCM group ##
################################################
# note there may be some internal bugs for the data labeled as NaN
# https://github.com/cBioPortal/cgdsr/issues/2
getmutations <- function() {
  mutations <- getProfileData(mycgds,
                              c("BRAF", "NRAS", "AKT", "TP53"),
                              "skcm_tcga_pan_can_atlas_2018_mutations",
                              "skcm_tcga_pan_can_atlas_2018_all")
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

# CNV.data <- getCNV("BRAF")
CNV.data <- getCNV(c("BRAF", "NRAS", "PTEN"))

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
##  plot OS curve with clinic and mutation data from SKCM group ##
##################################################################
plotOS <- function(ge) {
  mycancerstudy <- getCancerStudies(mycgds)[19, 1]
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

plotOS("BRAF.mutations")
mutation.list <- c("BRAF.mutations",
                   "NRAS.mutations",
                   "AKT1.mutations",
                   "TP53.mutations")
names(mutation.list) <- mutation.list
sapply(mutation.list, plotOS)

#########################################################################
##  plot OS curve with clinic and DNFA expression data from SKCM group ##
#########################################################################
## PANCAN dataset from cBioportal does not offer clinical results for OS
## need to the following script uses OS data from TCGA provisional data
## and combine OS data with the expression data from pancan study
plotDNFAOS <- function(DNFA) {
  mycancerstudy <- getCancerStudies(mycgds)[198, 1]  # "skcm_tcga_pan_can_atlas_2018"
  mycaselist <- getCaseLists(mycgds, mycancerstudy)[4, 1] # "All tumor samples (448 samples)"
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

plotDNFAOS("ACACA")
DNFA.list <- c("ACACA", "SCD", "ACLY", "FASN", "SREBF1", "MITF")
names(DNFA.list) <- DNFA.list
sapply(DNFA.list, plotDNFAOS)


##############################################################################
##  plot OS curve with clinic and DNFA expression data from all TCGA groups ##
##############################################################################
## PANCAN dataset from cBioportal does not offer clinical results for OS
## need to the following script uses OS data from TCGA provisional data
## and combine OS data with the expression data from pancan study
plotDNFAOS <- function(DNFA) {
  mycancerstudy <- getCancerStudies(mycgds)[201, 1]  # "skcm_tcga_pan_can_atlas_2018"
  mycaselist <- getCaseLists(mycgds, mycancerstudy)[4, 1] # "All tumor samples (448 samples)"
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

plotDNFAOS("ACACA")
DNFA.list <- c("ACACA", "SCD", "ACLY", "FASN", "SREBF1", "MITF")
names(DNFA.list) <- DNFA.list
sapply(DNFA.gene, plotDNFAOS)


##############################################################################
## Kaplan-Meier curve with clinic and DNFA RNAseq data from all TCGA groups ##
##############################################################################
plot.km.all.tcga <- function(DNFA) {
  #  mycgds <- CGDS("http://www.cbioportal.org/")
  #  test(mycgds)
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  pan.tcga.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, PanCancer Atlas)",
         getCancerStudies(mycgds)$name), ]
  pan.tcga.study.list <- pan.tcga.studies$cancer_study_id
  names(pan.tcga.study.list) <- pan.tcga.study.list
  pan.tcga.caselist <- lapply(pan.tcga.study.list,
                              caselist)
  pan.tcga.geneticprofile <- lapply(pan.tcga.study.list,
                                    geneticprofile)
  pan.caselist.RNAseq <- function(x) {
    pan.tcga.caselist[[x]][grep("tcga_pan_can_atlas_2018_rna_seq_v2_mrna", 
                                pan.tcga.caselist[[x]]$case_list_id), ][1, 1]
  } # pancancer group does not contain OS data
  pan.geneticprofile.RNAseq <- function(x) {
    pan.tcga.geneticprofile[[x]][
      grep("mRNA Expression, RSEM",
           pan.tcga.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
  }
  tcga.profiledata.RNAseq <- function(genename,
                                      geneticprofile,
                                      caselist) {
    getProfileData(mycgds,
                   genename,
                   geneticprofile,
                   caselist)
  }
  pan.tcga.gene.RNAseq <- function(x, y) {
    tcga.profiledata.RNAseq(x,
                            pan.geneticprofile.RNAseq(y),
                            pan.caselist.RNAseq(y))
  }
  pan.tcga.DNFA.RNAseq <- function(y) {
    pan.tcga.gene.RNAseq(x = DNFA, y)
  }
  ## test1 <- SCD.tcga.RNAseq(y = "skcm_tcga")
  ## test ## try to keep patient ID in the rowname
  test2 <- lapply(pan.tcga.study.list, pan.tcga.DNFA.RNAseq)
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
  pro.tcga.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)", getCancerStudies(mycgds)$name), ]
  # "pro.tcga.study.list" contains all the tcga provisional cancer studies
  pro.tcga.study.list <- pro.tcga.studies$cancer_study_id
  names(pro.tcga.study.list) <- pro.tcga.study.list
  # three datasets donot have OS data and cause bugs remove them
  bug.data.set <- names(pro.tcga.study.list) %in% c("meso_tcga", "pcpg_tcga", "ucs_tcga")
  pro.tcga.study.list <- pro.tcga.study.list[!bug.data.set]
  all.tcga.clinic.data <- lapply(pro.tcga.study.list, tcga.clinic.data)
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
    ggplot2::autoplot(km,
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

