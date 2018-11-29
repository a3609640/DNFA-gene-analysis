install.packages("cgdsr")
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
tcga.pro.studies <- getCancerStudies(mycgds)[grep("(TCGA, Provisional)", getCancerStudies(mycgds)$name), ]
# "tcag_study_list" is a vector containing all the tcga cancer studies that I would to analyze for DNFA gene expression
tcga.study.list <- tcga.pro.studies$cancer_study_id
names(tcga.study.list) <- tcga.study.list

caselist <- function(x) getCaseLists(mycgds, x)
geneticprofile <- function(x) getGeneticProfiles(mycgds, x)
# use lappy to pull out all the caselists within tcga.study.list
# because we named each elements in tcga.study.list (names(tcga.study.list) <- tcga.study.list),
# lappy will return a large list, each element (with a cancer study name) in that list is a data-table
tcga.pro.caselist <- lapply(tcga.study.list, 
                            caselist)
tcga.pro.geneticprofile <- lapply(tcga.study.list, 
                                  geneticprofile)
# for example, tcga.pro.caselist[[1]] shows the dataframe of caselist in laml study group.
# we want to choose case_list_id that is labeled with laml_tcga_rna_seq_v2_mrna, tcag_provisional_caselist[[1][8,1]
# a <- tcga.pro.caselist[[1]][grep("tcga_rna_seq_v2_mrna",
#                                       tcag_provisional_caselist[[1]]$case_list_id), ][1,1]
# b <- tcga.pro.geneticprofile[[1]][grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
# double backslash \\ suppress the special meaning of ( ) in regular expression
#                                             tcga.pro.geneticprofile[[1]]$genetic_profile_name), ][1,1]
# how do we do this for all study groups from [[1]] to  [[32]]?
caselist.RNAseq <- function(x) {
  tcga.pro.caselist[[x]][grep("tcga_rna_seq_v2_mrna", tcga.pro.caselist[[x]]$case_list_id), ][1, 1]
}

geneticprofile.RNAseq <- function(x) {
  tcga.pro.geneticprofile[[x]][grep("mRNA expression \\(RNA Seq V2 RSEM\\)", tcga.pro.geneticprofile[[x]]$genetic_profile_name), ][1, 1]
}
# test the functions: caselist.RNAseq () and geneticprofile.RNAseq ()
# caselist.RNAseq = caselist.RNAseq ('acc_tcga')
# geneticprofile.RNAseq = geneticprofile.RNAseq ('acc_tcga')
# We wrap two functions: geneticprofile.RNAseq(x), caselist.RNAseq(x) within TCGA_ProfileData_RNAseq(x)
tcga.profiledata.RNAseq <- function(genename,
                                    geneticprofile,
                                    caselist) {
  getProfileData(mycgds, genename, geneticprofile, caselist)
}

SCD.tcga.RNAseq <- function(x) {
  tcga.profiledata.RNAseq("SCD", geneticprofile.RNAseq(x), caselist.RNAseq(x))
}
# test the wrapping function, then use each element in tcga_study_list vector as the x in SCD_TCGA_RNAseq(x)
# data5 = TCGA_RNAseq ('acc_tcga')
SCD.tcga.RNAseq.all.studies <- lapply(tcga.study.list, SCD.tcga.RNAseq)

# use the melt function from reshape2 package.
df <- melt(SCD.tcga.RNAseq.all.studies)
# Separate boxplots for each data.frame
qplot(factor(L1), log2(value), data = df, geom = "boxplot")
# boxplot graph for SCD gene expression across all tumor types, and order the X axis based on the gene expression level
mean <- within(df, L1 <- reorder(L1, log2(value), median))
ggplot(mean, 
       aes(x    = L1, 
           y    = log2(value), 
           fill = L1)) +
  geom_boxplot(alpha    = .01, 
               width    = .5, 
               position = position_dodge(width = .9)) +
  theme_bw() +
  labs(x = "Tumor types (TCGA)", 
       y = paste("log2(SCD RNA counts)")) +
  theme(axis.title  = element_text(face   = "bold", 
                                   size   = 9, 
                                   color  = "black"), 
        axis.text   = element_text(size   = 9, 
                                   angle  = 45, 
                                   hjust  = 1, 
                                   face   = "bold", 
                                   color  = "black"), 
        axis.line.x = element_line(color  = "black"), 
        axis.line.y = element_line(color  = "black"), 
        panel.grid  = element_blank(), 
        strip.text  = element_text(face   = "bold", 
                                   size   = 9, 
                                   colour = "black"), 
        legend.position = "none")

##############################################
## Get DNFA gene expression from SKCM group ##
##############################################
skcm_case <- getCaseLists(mycgds, "skcm_tcga")
skcm_tcga_all <- getCaseLists(mycgds, "skcm_tcga")[2, 1]

# Get available genetic profiles
SKCMgeneticprofile <- getGeneticProfiles(mycgds, "skcm_tcga")
DNFA.RNAseq <- getProfileData(mycgds, 
                              c("ACACA", "FASN", "SCD", "ACLY", "ACSS2", "ACSL1", "LDLR", 
                                "SREBF1", "SREBF2","MITF", "BRAF", "NRAS", "PTEN"),
                              "skcm_tcga_rna_seq_v2_mrna",
                              "skcm_tcga_all")

################################################
## Get oncogene mutation data from SKCM group ##
################################################
getmutations <- function (x) {
  mutations <- getProfileData(mycgds, 
                            c("BRAF","NRAS","AKT","TP53"), 
                            "skcm_tcga_mutations", 
                            "skcm_tcga_all")
  colnames(mutations) <- paste0(colnames(mutations), '.mutations')
  v <- rownames(mutations)
  # each mutation column contains three types of data: 
  # mutation (V600E), NAN (wildtype), NA (not sequenced).
  relabel.mutations <- function(gene) {
    mutations[, gene] <- ifelse(mutations[, gene] == "NaN", "Wildtype", "Mutated")
    }
# use sapply , input as a matrix, and output as a matrix too.
  mutations <- sapply(colnames(mutations), relabel.mutations)
  mutations <- as.data.frame(mutations)
# the sapply function return a new matrix lacking row names. need to add that back
  mutations2 <- cbind(Row.Names = v, mutations)
# mutations <- mutations2[ , -1]
  rownames(mutations) <- mutations2[ , 1]
  return (mutations)
  }
mutations <- getmutations (c("BRAF", "NRAS", "AKT", "TP53"))
attributes(mutations)

###########################################
## Get oncogene CNV data from SKCM group ##
###########################################
getCNV <- function (x) {
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
  # use sapply , input as a matrix, and output as a matrix too.
  CNV <- sapply(colnames(CNV), factor.CNV)
  CNV <- as.data.frame(CNV)
  CNV2 <- cbind(Row.Names = v, CNV)
  # mutations <- mutations2[ , -1]
  rownames(CNV) <- CNV2[ , 1]
  return (CNV)
}
CNV <- getCNV (c("BRAF", "NRAS", "PTEN"))
attributes(CNV)

###############################################################################
## examine the correlation between mutation status and DNFA expression level ##
###############################################################################
mutations.DNFA.RNAseq <- cbind(mutations, DNFA.RNAseq)
mutations.DNFA.RNAseq <- na.omit(mutations.DNFA.RNAseq)

plot.mutations.RNAseq <- function (genemutations) {
  ggplot(mutations.DNFA.RNAseq, 
         aes(x     = mutations.DNFA.RNAseq[ ,genemutations], # refer variable in a function
             y     = log2(SCD), 
             color = mutations.DNFA.RNAseq[ ,genemutations])) + 
    geom_boxplot(alpha = .01, 
                 width = .5) +
    labs(x = genemutations, 
         y = paste("log2(", "SCD", "RNA counts)"))+ 
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
          legend.position = "none")+
    geom_dotplot(binaxis  = "y", # stack the dots along the y-axis and group them along x-axis
               binwidth   = .1,
               stackdir   = "center",
               fill       = NA)
  }

plot.mutations.RNAseq("NRAS.mutations")
lapply(colnames(mutations), plot.mutations.RNAseq)

########################################################################################
##  check the distribution of RNA-seq, then we choose the stastics comparison methods ##
########################################################################################
BRAF.Mutated <- mutations.DNFA.RNAseq$SCD[mutations.DNFA.RNAseq$BRAF.mutations == "Mutated"]
BRAF.Wildtype <- mutations.DNFA.RNAseq$SCD[mutations.DNFA.RNAseq$BRAF.mutations == "Wildtype"]
plot(density(BRAF.Mutated))
plot(density(BRAF.Wildtype))
# Normality test
shapiro.test(BRAF.Mutated) # p-value < 2.2e-16
shapiro.test(BRAF.Wildtype) # p-value < 2.2e-16
# Welch two sample t-test
t.test(BRAF.Mutated, BRAF.Wildtype) # p-value = 0.5361

NRAS.Mutated <- mutations.DNFA.RNAseq$SCD[mutations.DNFA.RNAseq$NRAS.mutations == "Mutated"]
NRAS.Wildtype <- mutations.DNFA.RNAseq$SCD[mutations.DNFA.RNAseq$NRAS.mutations == "Wildtype"]
plot(density(NRAS.Mutated))
plot(density(NRAS.Wildtype))
# Normality test
shapiro.test(NRAS.Mutated) # p-value = 2.279e-13
shapiro.test(NRAS.Wildtype) # p-value < 2.2e-16
t.test(NRAS.Mutated, NRAS.Wildtype) # p-value = 0.9444

AKT.Mutated <- mutations.DNFA.RNAseq$SCD[mutations.DNFA.RNAseq$AKT.mutations == "Mutated"]
AKT.Wildtype <- mutations.DNFA.RNAseq$SCD[mutations.DNFA.RNAseq$AKT.mutations == "Wildtype"]
plot(density(AKT.Mutated))
plot(density(AKT.Wildtype))
# Normality test
shapiro.test(AKT.Mutated) # p-value = 0.8429 not a normal distribution
shapiro.test(AKT.Wildtype) # p-value < 2.2e-16
# t.test(AKT.Mutated, AKT.Wildtype)
# use non-parametric test
wilcox.test(AKT.Mutated, AKT.Wildtype) # p-value = 0.434

TP53.Mutated <- mutations.DNFA.RNAseq$SCD[mutations.DNFA.RNAseq$TP53.mutations == "Mutated"]
TP53.Wildtype <- mutations.DNFA.RNAseq$SCD[mutations.DNFA.RNAseq$TP53.mutations == "Wildtype"]
plot(density(TP53.Mutated))
plot(density(TP53.Wildtype))
# Normality test
shapiro.test(TP53.Mutated) # p-value = 0.009181
shapiro.test(TP53.Wildtype) # p-value < 2.2e-16
# Welch Two Sample t-test
t.test(TP53.Mutated, TP53.Wildtype) # p-value = 0.06594
# Wilcoxon rank sum test with continuity correction
wilcox.test(TP53.Mutated, TP53.Wildtype) # p-value = 0.761

#############################################################################
## compare the correlation between oncogene CNV and DNFA gene RNA-seq data ##
#############################################################################
CNV.DNFA.RNAseq <- cbind(CNV, DNFA.RNAseq)
# na.omit will eleminate NaN as well
CNV.DNFA.RNAseq <- na.omit(CNV.DNFA.RNAseq)
# make a large function to plot all genes
plot.CNV.RNAseq <- function (geneCNV, RNAseq) {
  print(ggplot(CNV.DNFA.RNAseq, 
               aes(x     = CNV.DNFA.RNAseq[ ,geneCNV], 
                   y     = log2(CNV.DNFA.RNAseq[ ,RNAseq]), 
                   color = CNV.DNFA.RNAseq[ ,geneCNV])) + 
          geom_boxplot(alpha      = .01, 
                       width      = .5) + 
          geom_dotplot(binaxis    = "y", 
                       binwidth   = .1, 
                       stackdir   = "center", 
                       fill       = NA) + 
          scale_x_discrete(limits = c("-2", "-1", "0", "1", "2"), # do not show NaN data
                           labels = c("Homdel", "Hetlos", "Diploid", "Gain", "Amp"), 
                           drop   = FALSE) + 
          labs(x = geneCNV, 
               y = paste("log2(", RNAseq, "RNA counts)"))+ 
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
                legend.position = "none"))
  }

plot.CNV.RNAseq ("NRAS.CNV", "BRAF")
plot.CNV.RNAseq ("NRAS.CNV", "SCD")
# however, mapply function only works with print function in ggplot.
# three combinations: BRAF.CNV, BRAF; NRAS.CNV, NRAS; PTEN.CNV, PTEN
mapply(plot.CNV.RNAseq, 
       geneCNV = c("BRAF.CNV","NRAS.CNV", "PTEN.CNV"), 
       RNAseq = c("BRAF", "NRAS", "PTEN"))
# use all combinations of geneCNV and RNAseq
sapply(c("BRAF.CNV","NRAS.CNV", "PTEN.CNV"), function(x) mapply(plot.CNV.RNAseq,x,c("BRAF", "NRAS", "PTEN", "SCD")))

# compare the variance among different groups using a nonparametri test
fligner.test(SCD ~ BRAF.CNV, data = BRAF.CNV.DNFA.RNAseq) # p-value = 0.179
fligner.test(BRAF ~ BRAF.CNV, data = BRAF.CNV.DNFA.RNAseq) # p-value = 0.02785
#################################################################################################

# compare the variance among different groups using a nonparametri test
fligner.test(SCD ~ NRAS.CNV, data = NRAS.CNV.DNFA.RNAseq) # p-value = 0.1973
fligner.test(NRAS ~ NRAS.CNV, data = NRAS.CNV.DNFA.RNAseq) # p-value = 1.9e-08
#################################################################################################
# compare the variance among different groups using a nonparametri test
fligner.test(SCD ~ PTEN.CNV, data = PTEN.CNV.DNFA.RNAseq) # p-value = 0.0004647 (pten loss correlates with scd decrease?)
fligner.test(PTEN ~ PTEN.CNV, data = PTEN.CNV.DNFA.RNAseq) # p-value = p-value = 5.948e-08

# statistic comparision of SCD expression between PTEN hetero deletion and PTEN diploid
PTEN.Hetlos <- PTEN.CNV.DNFA.RNAseq$SCD[PTEN.CNV.DNFA.RNAseq$PTEN.CNV == "-1"]
PTEN.Diploid <- PTEN.CNV.DNFA.RNAseq$SCD[PTEN.CNV.DNFA.RNAseq$PTEN.CNV == "0"]
plot(density(PTEN.Hetlos))
plot(density(PTEN.Diploid))
# Normality test
shapiro.test(PTEN.Hetlos) # p-value = 6.252e-08
shapiro.test(PTEN.Diploid) # p-value = 5.124e-13
# Welch Two Sample t-test
t.test(PTEN.Hetlos, PTEN.Diploid) # p-value = 0.0008492 (pten loss correlates with scd decrease)
####################################################################################################################
####################################################################################################################


mycancerstudy <- getCancerStudies(mycgds)[193, 1]
mycaselist <- getCaseLists(mycgds, mycancerstudy)[1, 1]
skcm.clinicaldata <- getClinicalData(mycgds, 
                                     mycaselist)

skcm.clinicaldata$rn <- rownames(skcm.clinicaldata)
BRAF.mutations.DNFA.RNAseq$rn <- rownames(BRAF.mutations.DNFA.RNAseq)
NRAS.mutations.DNFA.RNAseq$rn <- rownames(NRAS.mutations.DNFA.RNAseq)
TP53.mutations.DNFA.RNAseq$rn <- rownames(TP53.mutations.DNFA.RNAseq)
BRAF.CNV.DNFA.RNAseq$rn <- rownames(BRAF.CNV.DNFA.RNAseq)
NRAS.CNV.DNFA.RNAseq$rn <- rownames(NRAS.CNV.DNFA.RNAseq)
PTEN.CNV.DNFA.RNAseq$rn <- rownames(PTEN.CNV.DNFA.RNAseq)

df <- join_all(list(skcm.clinicaldata[c("OS_MONTHS", 
                                        "OS_STATUS", 
                                        'rn')],
                    BRAF.mutations.DNFA.RNAseq[c("BRAF.mutations", 
                                                 'rn')],
                    NRAS.mutations.DNFA.RNAseq[c("NRAS.mutations", 
                                                 'rn')],
                    TP53.mutations.DNFA.RNAseq[c("TP53.mutations", 
                                                 'rn')],
                    BRAF.CNV.DNFA.RNAseq[c("BRAF.CNV", 
                                           'rn')],
                    NRAS.CNV.DNFA.RNAseq[c("NRAS.CNV", 
                                           'rn')],
                    PTEN.CNV.DNFA.RNAseq[c("PTEN.CNV", 
                                           'rn')]), 
               by   = 'rn', 
               type = 'full')

df$SurvObj <- with(df, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))
km.by.mutation <- survfit(SurvObj ~ NRAS.mutations, 
                              data      = df, 
                              conf.type = "log-log")
autoplot(km.by.mutation)
## Load survival package


