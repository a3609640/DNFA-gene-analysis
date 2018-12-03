# install.packages("cgdsr")
library(car)
library(cgdsr)
library(ggfortify)
library(ggplot2)
library(plyr)
library(reshape2)
library(survival)

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
  tcga.pro.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)", getCancerStudies(mycgds)$name), ]
  # "tcag_study_list" is a vector containing all the tcga cancer studies 
  # that I would to analyze for DNFA gene expression
  tcga.study.list <- tcga.pro.studies$cancer_study_id
  names(tcga.study.list) <- tcga.study.list
  caselist <- function(x) getCaseLists(mycgds, x)
  geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
  # use lappy to pull out all the caselists within tcga.study.list
  # because we named each elements in tcga.study.list 
  # (names(tcga.study.list) <- tcga.study.list),
  # lappy will return a large list, each element (with a cancer study name) 
  # in that list is a data-table
  tcga.pro.caselist <- lapply(tcga.study.list, caselist)
  tcga.pro.geneticprofile <- lapply(tcga.study.list, geneticprofile)
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
    tcga.pro.caselist[[x]][
      grep("tcga_rna_seq_v2_mrna", tcga.pro.caselist[[x]]$case_list_id), ][1, 1]
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
lapply(DNFA.gene, plotDNFA)

##############################################
## Get DNFA gene expression from SKCM group ##
##############################################
# skcm_case <- getCaseLists(mycgds, "skcm_tcga")
# skcm_tcga_all <- getCaseLists(mycgds, "skcm_tcga")[2, 1]
# Get available genetic profiles
# SKCMgeneticprofile <- getGeneticProfiles(mycgds, "skcm_tcga")
DNFA.RNAseq.data <- getProfileData(mycgds,
                                   c("ACACA", "FASN", "SCD", "ACLY", "ACSS2",
                                     "ACSL1", "LDLR", "SREBF1", "SREBF2",
                                     "MITF", "BRAF", "NRAS", "PTEN"),
                                   "skcm_tcga_rna_seq_v2_mrna",
                                   "skcm_tcga_all")

################################################
## Get oncogene mutation data from SKCM group ##
################################################
getmutations <- function(x) {
  mutations <- getProfileData(mycgds,
                              c("BRAF","NRAS","AKT","TP53"),
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
  print(attributes(mutations))
  return(mutations)
  }
mutations.data <- getmutations(c("BRAF", "NRAS", "AKT", "TP53"))

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
  print(attributes(CNV))
  return(CNV)
}

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
    geom_dotplot(binaxis    = "y", 
#stack the dots along the y-axis and group them along x-axis                 
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
# statistic comparision of SCD expression between PTEN heterdeletion and diploid
# (?pten loss correlates with scd decrease?)

############################################
##  Retrieve clinic data from SKCM group  ##
############################################
mycancerstudy <- getCancerStudies(mycgds)[193, 1]
mycaselist <- getCaseLists(mycgds, mycancerstudy)[1, 1]
skcm.clinicaldata <- getClinicalData(mycgds, mycaselist)
skcm.clinicaldata$rn <- rownames(skcm.clinicaldata)

mutations.DNFA.RNAseq <- cbind(mutations.data, DNFA.RNAseq.data)
mutations.DNFA.RNAseq <- na.omit(mutations.DNFA.RNAseq)
mutations.DNFA.RNAseq$rn <- rownames(mutations.DNFA.RNAseq)

CNV.DNFA.RNAseq <- cbind(CNV.data, DNFA.RNAseq.data)
# na.omit cannot eleminate NaN here!
toBeRemoved <- which(CNV.DNFA.RNAseq$BRAF.CNV == "NaN")
CNV.DNFA.RNAseq <- CNV.DNFA.RNAseq[-toBeRemoved,]
CNV.DNFA.RNAseq$rn <- rownames(CNV.DNFA.RNAseq)

df <- join_all(list(skcm.clinicaldata[c("OS_MONTHS", 
                                        "OS_STATUS", 
                                        'rn')],
                    mutations.DNFA.RNAseq[c("BRAF.mutations", 
                                            "NRAS.mutations",
                                            "TP53.mutations",
                                            'rn')],
                    CNV.DNFA.RNAseq[c("BRAF.CNV",
                                      "NRAS.CNV",
                                      "PTEN.CNV",
                                      'rn')]),
               by   = 'rn', 
               type = 'full')

df$SurvObj <- with(df, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))
km.by.mutation <- survfit(SurvObj ~ NRAS.mutations,
                          data      = df,
                          conf.type = "log-log")
autoplot(km.by.mutation)
library("survminer")
ggsurvplot(km.by.mutation, data = df)
