install.packages('cgdsr')
library(cgdsr)
library(ggplot2)
library(reshape2)

# Create CGDS object
mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
test(mycgds)
# Get list of cancer studies at server
getCancerStudies(mycgds)
# Get cases from TCGA provisional studies only
tcga_provisional_studies <- getCancerStudies(mycgds)[grep("(TCGA, Provisional)", 
                                                         getCancerStudies(mycgds)$name), ]
# "tcag_study_list" is a vector containing all the tcga cancer studies that I would to analyze for DNFA gene expression
tcga_study_list <- tcga_provisional_studies$cancer_study_id
names(tcga_study_list) <- tcga_study_list

caselist <- function (x) getCaseLists(mycgds, x)
geneticprofile <- function (x) getGeneticProfiles(mycgds, x)
# use lappy to pull out all the caselists within tcag_study_list
# because we named each elements in tcga_study_list (names(tcag_study_list) <- tcag_study_list), 
# lappy will return a large list, each element (with a cancer study name) in that list is a data-table
tcga_provisional_caselist <- lapply (tcga_study_list, caselist)
tcga_provisional_geneticprofile <- lapply (tcga_study_list, geneticprofile)
# for example, tcag_provisional_caselist[[1]] shows the dataframe of caselist in laml study group.
# we want to choose case_list_id that is labeled with laml_tcga_rna_seq_v2_mrna, tcag_provisional_caselist[[1][8,1]
  # a =tcag_provisional_caselist[[1]][grep("tcga_rna_seq_v2_mrna",
  #                                       tcag_provisional_caselist[[1]]$case_list_id), ][1,1]
  # b =tcag_provisional_geneticprofile[[1]][grep("mRNA expression \\(RNA Seq V2 RSEM\\)", 
  # double backslash \\ suppress the special meaning of ( ) in regular expression
  #                                             tcag_provisional_geneticprofile[[1]]$genetic_profile_name), ][1,1]
# how do we do this for all study groups from [[1]] to  [[32]]?
caselist_RNAseq <- function(x) {
                                tcga_provisional_caselist[[x]][grep("tcga_rna_seq_v2_mrna", 
                                                                    tcga_provisional_caselist[[x]]$case_list_id),
                                                               ][1,1]
                                }

geneticprofile_RNAseq <- function(x) {
                                      tcga_provisional_geneticprofile[[x]][grep("mRNA expression \\(RNA Seq V2 RSEM\\)", 
                                                                                tcga_provisional_geneticprofile[[x]]$genetic_profile_name),
                                                                           ][1,1]
                                      }
# test the functions: caselist_RNAseq () and geneticprofile_RNAseq ()
    # caselist_RNAseq = caselist_RNAseq ('acc_tcga')
    # geneticprofile_RNAseq = geneticprofile_RNAseq ('acc_tcga')
# We wrap two functions: geneticprofile_RNAseq(x), caselist_RNAseq(x) within TCGA_ProfileData_RNAseq(x)
TCGA_ProfileData_RNAseq <- function(genename, geneticprofile, caselist) {
                                                                         getProfileData(mycgds, 
                                                                                        genename, 
                                                                                        geneticprofile, 
                                                                                        caselist)
                                                                         }

SCD_TCGA_RNAseq <- function(x) {
                                TCGA_ProfileData_RNAseq('SCD', 
                                geneticprofile_RNAseq(x), 
                                caselist_RNAseq(x))
                                }  
# test the wrapping function, then use each element in tcga_study_list vector as the x in SCD_TCGA_RNAseq(x)
  # data5 = TCGA_RNAseq ('acc_tcga')
SCD_TCGA_RNAseq_all_Studies <- lapply (tcga_study_list, SCD_TCGA_RNAseq)

# use the melt function from reshape2 package. 
df <- melt(SCD_TCGA_RNAseq_all_Studies)
# Separate boxplots for each data.frame
qplot(factor(L1), log2(value), data = df, geom = "boxplot")
# boxplot graph for SCD gene expression across all tumor types, and order the X axis based on the gene expression level
mean <- within(df, L1 <- reorder(L1, log2(value), median))
ggplot(mean, aes(x    = L1, 
                 y    = log2(value), 
                 fill = L1)) + 
  geom_boxplot(alpha    = .01, 
               width    = .5, 
               position = position_dodge(width = .9)) +
  theme_bw() +
  labs(x = "Tumor types (TCGA)",
       y = paste("log2(SCD RNA counts)")) +
  theme(axis.title      = element_text(face = "bold",size = ,color =" black"),
        axis.text       = element_text(size = 9, angle = 45, hjust = 1, face = "bold", color = "black"),
        axis.line.x     = element_line(color = "black"),
        axis.line.y     = element_line(color = "black"),
        panel.grid      = element_blank(),
        strip.text      = element_text(face = "bold", size = 9, colour = "black"),
        legend.position = "none")

#################################################################################################
## get DNFA gene expression from SKCM group and check the correlation with oncogenic mutations ##
#################################################################################################
skcm_case <- getCaseLists(mycgds,'skcm_tcga')
skcm_tcga_all <- getCaseLists(mycgds,'skcm_tcga')[2,1]

# Get available genetic profiles
SKCMgeneticprofile <- getGeneticProfiles(mycgds,'skcm_tcga')
DNFA.RNAseq=getProfileData(mycgds,
                           c('ACACA','FASN','SCD','ACLY','ACSS2','ACSL1','LDLR','SREBF1','SREBF2',
                             'MITF','BRAF', 'NRAS', 'PTEN'),
                           'skcm_tcga_rna_seq_v2_mrna',
                           'skcm_tcga_all')

#################################################################################################
# Get data slices for a specified list of genes, genetic profile and case list
BRAF.mutations=getProfileData(mycgds, 'BRAF',
                              'skcm_tcga_mutations',
                              'skcm_tcga_all')
colnames(BRAF.mutations) <- c("BRAF.mutations")
BRAF.mutations.DNFA.RNAseq <- cbind(BRAF.mutations, DNFA.RNAseq)
levels(BRAF.mutations.DNFA.RNAseq$BRAF.mutations) <- c(levels(BRAF.mutations.DNFA.RNAseq$BRAF.mutations), 
                                                       "Mutated", 
                                                       "Wildtype")
BRAF.mutations.DNFA.RNAseq$BRAF.mutations[BRAF.mutations.DNFA.RNAseq$BRAF.mutations!='NA'] <- "Mutated"
BRAF.mutations.DNFA.RNAseq$BRAF.mutations[is.na(BRAF.mutations.DNFA.RNAseq$BRAF.mutations)] <- "Wildtype" 
BRAF.mutations.DNFA.RNAseq <- na.omit(BRAF.mutations.DNFA.RNAseq)
ggplot(BRAF.mutations.DNFA.RNAseq, 
       aes(x     = BRAF.mutations, 
           y     = log2(SCD), 
           color = BRAF.mutations)) + 
  geom_boxplot(alpha = .01, 
               width = .5) +
  geom_dotplot(binaxis  = "y", # stack the dots along the y-axis and group them along x-axis
               binwidth = .1, 
               stackdir = "center", 
               fill     = NA)

# density plot to first check the distribution of RNA-seq, then we choose the stastics comparison method
BRAF.Mutated = BRAF.mutations.DNFA.RNAseq$SCD[BRAF.mutations.DNFA.RNAseq$BRAF.mutations == "Mutated"]
BRAF.Wildtype = BRAF.mutations.DNFA.RNAseq$SCD[BRAF.mutations.DNFA.RNAseq$BRAF.mutations == "Wildtype"]
plot(density(BRAF.Mutated))
plot(density(BRAF.Wildtype))
# Normality test
shapiro.test(BRAF.Mutated) # p-value < 2.2e-16
shapiro.test(BRAF.Wildtype) # p-value < 2.2e-16
# Welch two sample t-test
t.test(BRAF.Mutated, BRAF.Wildtype) # p-value = 0.5361

#################################################################################################
NRAS.mutations=getProfileData(mycgds, 'NRAS',
                              'skcm_tcga_mutations',
                              'skcm_tcga_all')
colnames(NRAS.mutations) <- c("NRAS.mutations")
NRAS.mutations.DNFA.RNAseq <- cbind(NRAS.mutations, DNFA.RNAseq)
levels(NRAS.mutations.DNFA.RNAseq$NRAS.mutations) <- c(levels(NRAS.mutations.DNFA.RNAseq$NRAS.mutations), 
                                                       "Mutated", 
                                                       "Wildtype")
NRAS.mutations.DNFA.RNAseq$NRAS.mutations[NRAS.mutations.DNFA.RNAseq$NRAS.mutations!='NA'] <- "Mutated"
NRAS.mutations.DNFA.RNAseq$NRAS.mutations[is.na(NRAS.mutations.DNFA.RNAseq$NRAS.mutations)] <- "Wildtype" 
NRAS.mutations.DNFA.RNAseq = na.omit(NRAS.mutations.DNFA.RNAseq)
ggplot(NRAS.mutations.DNFA.RNAseq, 
       aes(x     = NRAS.mutations, 
           y     = log2(SCD), 
           color = NRAS.mutations)) + 
  geom_boxplot(alpha = .01, 
               width = .5) +
  geom_dotplot(binaxis  = "y", 
               binwidth = .1, 
               stackdir = "center", 
               fill     = NA)

# density plot to first check the distribution of RNA-seq, then we choose the stastics comparison method
NRAS.Mutated <- NRAS.mutations.DNFA.RNAseq$SCD[NRAS.mutations.DNFA.RNAseq$NRAS.mutations == "Mutated"]
NRAS.Wildtype <- NRAS.mutations.DNFA.RNAseq$SCD[NRAS.mutations.DNFA.RNAseq$NRAS.mutations == "Wildtype"]
plot(density(NRAS.Mutated))
plot(density(NRAS.Wildtype))
# Normality test
shapiro.test(NRAS.Mutated) # p-value = 2.279e-13
shapiro.test(NRAS.Wildtype) # p-value < 2.2e-16
t.test(NRAS.Mutated, NRAS.Wildtype) # p-value = 0.9444

#################################################################################################
AKT.mutations=getProfileData(mycgds, 
                             'AKT',
                             'skcm_tcga_mutations',
                             'skcm_tcga_all')
colnames(AKT.mutations) <- c("AKT.mutations")
AKT.mutations.DNFA.RNAseq=cbind(AKT.mutations, DNFA.RNAseq)
levels(AKT.mutations.DNFA.RNAseq$AKT.mutations) <- c(levels(AKT.mutations.DNFA.RNAseq$AKT.mutations), 
                                                     "Mutated", 
                                                     "Wildtype")
AKT.mutations.DNFA.RNAseq$AKT.mutations[AKT.mutations.DNFA.RNAseq$AKT.mutations!='NA'] <- "Mutated"
AKT.mutations.DNFA.RNAseq$AKT.mutations[is.na(AKT.mutations.DNFA.RNAseq$AKT.mutations)] <- "Wildtype" 
AKT.mutations.DNFA.RNAseq = na.omit(AKT.mutations.DNFA.RNAseq)
ggplot(AKT.mutations.DNFA.RNAseq, 
       aes(x     = AKT.mutations, 
           y     = log2(SCD), 
           color = AKT.mutations)) + 
  geom_boxplot(alpha = .01, 
               width = .5) +
  geom_dotplot(binaxis  = "y", 
               binwidth = .1, 
               stackdir = "center", 
               fill     = NA)

# density plot to first check the distribution of RNA-seq, then we choose the stastics comparison method
AKT.Mutated <- AKT.mutations.DNFA.RNAseq$SCD[AKT.mutations.DNFA.RNAseq$AKT.mutations == "Mutated"]
AKT.Wildtype <- AKT.mutations.DNFA.RNAseq$SCD[AKT.mutations.DNFA.RNAseq$AKT.mutations == "Wildtype"]
plot(density(AKT.Mutated)) 
plot(density(AKT.Wildtype))
# Normality test
shapiro.test(AKT.Mutated) # p-value = 0.8429 not a normal distribution
shapiro.test(AKT.Wildtype)
# t.test(AKT.Mutated, AKT.Wildtype)
# use non-parametric test
wilcox.test(AKT.Mutated, AKT.Wildtype) # p-value = 0.434

#################################################################################################
TP53.mutations=getProfileData(mycgds, 
                              'TP53',
                              'skcm_tcga_mutations',
                              'skcm_tcga_all')
colnames(TP53.mutations) <- c("TP53.mutations")
TP53.mutations.DNFA.RNAseq <- cbind(TP53.mutations, DNFA.RNAseq)
levels(TP53.mutations.DNFA.RNAseq$TP53.mutations) <- c(levels(TP53.mutations.DNFA.RNAseq$TP53.mutations), 
                                                       "Mutated", 
                                                       "Wildtype")
TP53.mutations.DNFA.RNAseq$TP53.mutations[TP53.mutations.DNFA.RNAseq$TP53.mutations!='NA'] <- "Mutated"
TP53.mutations.DNFA.RNAseq$TP53.mutations[is.na(TP53.mutations.DNFA.RNAseq$TP53.mutations)] <- "Wildtype" 
TP53.mutations.DNFA.RNAseq <- na.omit(TP53.mutations.DNFA.RNAseq)
ggplot(TP53.mutations.DNFA.RNAseq, 
       aes(x     = TP53.mutations, 
           y     = log2(SCD), 
           color = TP53.mutations)) + 
  geom_boxplot(alpha = .01, 
               width = .5) +
  geom_dotplot(binaxis  = "y", 
               binwidth = .1, 
               stackdir = "center", 
               fill     = NA)

# density plot to first check the distribution of RNA-seq, then we choose the stastics comparison method
TP53.Mutated <- TP53.mutations.DNFA.RNAseq$SCD[TP53.mutations.DNFA.RNAseq$TP53.mutations == "Mutated"]
TP53.Wildtype <- TP53.mutations.DNFA.RNAseq$SCD[TP53.mutations.DNFA.RNAseq$TP53.mutations == "Wildtype"]
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
BRAF.CNV <- getProfileData(mycgds, 
                           'BRAF',
                           "skcm_tcga_gistic",
                           "skcm_tcga_all")
colnames(BRAF.CNV) <- c("BRAF.CNV")
BRAF.CNV.DNFA.RNAseq <- cbind(BRAF.CNV, DNFA.RNAseq)
# there are three different ways to remove NaN in the dataset
BRAF.CNV.DNFA.RNAseq <- BRAF.CNV.DNFA.RNAseq[!is.nan(BRAF.CNV.DNFA.RNAseq$BRAF.CNV), ]
  # BRAF.CNV.DNFA.RNAseq = BRAF.CNV.DNFA.RNAseq[complete.cases(BRAF.CNV.DNFA.RNAseq), ]
  # BRAF.CNV.DNFA.RNAseq = na.omit(BRAF.CNV.DNFA.RNAseq)
# change $BRAF column from numeric to factor
BRAF.CNV.DNFA.RNAseq$BRAF.CNV <- as.factor(BRAF.CNV.DNFA.RNAseq$BRAF.CNV)
# some data set might miss one or two CNV categery, but we want to present all  
levels(BRAF.CNV.DNFA.RNAseq$BRAF.CNV) <- c( "2", "1",  "0",  "-1", "-2")
# plot SCD RNAseq ~ BRAF CNV
ggplot(BRAF.CNV.DNFA.RNAseq, 
       aes(x     = BRAF.CNV, 
           y     = log2(SCD), 
           color = BRAF.CNV)) + 
  geom_boxplot(alpha = .01, 
               width = .5) +
  geom_dotplot(binaxis  = "y", 
               binwidth = .1, 
               stackdir = "center", 
               fill     = NA) +
  scale_x_discrete(breaks = c("-2", "-1",  "0",  "1", "2"), 
                   labels = c("Homdel", "Hetlos", "Diploid", "Gain", "Amp"), 
                   drop   = FALSE) +
  theme(legend.position   = "none")

# plot BRAF RNAseq ~ BRAF CNV: positive control
ggplot(BRAF.CNV.DNFA.RNAseq, 
       aes(x     = BRAF.CNV, 
           y     = log2(BRAF), 
           color = BRAF.CNV)) + 
  geom_boxplot(alpha = .01, 
               width = .5) +
  geom_dotplot(binaxis  = "y", 
               binwidth = .1, 
               stackdir = "center", 
               fill     = NA) +
  scale_x_discrete(breaks = c("-2", "-1",  "0",  "1", "2"), 
                   labels = c("Homdel", "Hetlos", "Diploid", "Gain", "Amp"), 
                   drop   = FALSE) +
  theme(legend.position   = "none")

# compare the variance among different groups using a nonparametri test
fligner.test(SCD~BRAF.CNV, data=BRAF.CNV.DNFA.RNAseq) # p-value = 0.179
fligner.test(BRAF~BRAF.CNV, data=BRAF.CNV.DNFA.RNAseq) # p-value = 0.02785
#################################################################################################
NRAS.CNV=getProfileData(mycgds, 
                        "NRAS",
                        "skcm_tcga_gistic",
                        "skcm_tcga_all")
colnames(NRAS.CNV) <- c("NRAS.CNV")
NRAS.CNV.DNFA.RNAseq <- cbind(NRAS.CNV, DNFA.RNAseq)
NRAS.CNV.DNFA.RNAseq <- NRAS.CNV.DNFA.RNAseq[!is.nan(NRAS.CNV.DNFA.RNAseq$NRAS.CNV), ]
NRAS.CNV.DNFA.RNAseq$NRAS.CNV <- as.factor(NRAS.CNV.DNFA.RNAseq$NRAS.CNV)
levels(NRAS.CNV.DNFA.RNAseq$NRAS.CNV) <- c( "-2", "-1",  "0",  "1", "2")

# plot SCD RNAseq ~ NRAS CNV
ggplot(NRAS.CNV.DNFA.RNAseq, 
       aes(x     = NRAS.CNV, 
           y     = log2(SCD), 
           color = NRAS.CNV)) + 
  geom_boxplot(alpha = .01, 
               width = .5) +
  geom_dotplot(binaxis  = "y", 
               binwidth = .1, 
               stackdir = "center", 
               fill     = NA) +
  scale_x_discrete(breaks = c("-2", "-1",  "0",  "1", "2"), 
                   labels = c("Homdel", "Hetlos", "Diploid", "Gain", "Amp"), 
                   drop   = FALSE) +
  theme(legend.position   = "none")

# plot NRAS RNAseq ~ NRAS CNV: positive control
ggplot(NRAS.CNV.DNFA.RNAseq, 
       aes(x     = NRAS.CNV, 
           y     = log2(NRAS), 
           color = NRAS.CNV)) + 
  geom_boxplot(alpha = .01, 
               width = .5) +
  geom_dotplot(binaxis  = "y", 
               binwidth = .1, 
               stackdir = "center", 
               fill     = NA) +
  scale_x_discrete(breaks = c("-2", "-1",  "0",  "1", "2"), 
                   labels = c("Homdel", "Hetlos", "Diploid", "Gain", "Amp"), 
                   drop   = FALSE) +
  theme(legend.position   = "none")

# compare the variance among different groups using a nonparametri test
fligner.test(SCD~NRAS.CNV, data=NRAS.CNV.DNFA.RNAseq) # p-value = 0.1973
fligner.test(NRAS~NRAS.CNV, data=NRAS.CNV.DNFA.RNAseq) # p-value = 1.9e-08
#################################################################################################
PTEN.CNV=getProfileData(mycgds, 
                        "PTEN",
                        "skcm_tcga_gistic",
                        "skcm_tcga_all")
colnames(PTEN.CNV) <- c("PTEN.CNV")
PTEN.CNV.DNFA.RNAseq <- cbind(PTEN.CNV, DNFA.RNAseq)
PTEN.CNV.DNFA.RNAseq <- PTEN.CNV.DNFA.RNAseq[!is.nan(PTEN.CNV.DNFA.RNAseq$PTEN.CNV), ]
PTEN.CNV.DNFA.RNAseq$PTEN.CNV <- as.factor(PTEN.CNV.DNFA.RNAseq$PTEN.CNV)
levels(PTEN.CNV.DNFA.RNAseq$PTEN.CNV) <- c( "-2", "-1", "0", "1", "2")
ggplot(PTEN.CNV.DNFA.RNAseq, 
       aes(x     = PTEN.CNV, 
           y     = log2(SCD), 
           color = PTEN.CNV)) + 
  geom_boxplot(alpha = .01, 
               width = .5) +
  geom_dotplot(binaxis  = "y", 
               binwidth = .1, 
               stackdir = "center", 
               fill     = NA) +
  scale_x_discrete(breaks = c("-2", "-1", "0", "1", "2"), 
                   labels = c("Homdel", "Hetlos", "Diploid", "Gain", "Amp"), 
                   drop   = FALSE) +
  theme(legend.position   = "none")

# plot PTEN RNAseq ~ PTEN CNV: positive control
ggplot(PTEN.CNV.DNFA.RNAseq, 
       aes(x     = PTEN.CNV, 
           y     = log2(PTEN), 
           color = PTEN.CNV)) + 
  geom_boxplot(alpha = .01, 
               width = .5) +
  geom_dotplot(binaxis  = "y", 
               binwidth = .1, 
               stackdir = "center", 
               fill     = NA) +
  scale_x_discrete(breaks = c("-2", "-1", "0", "1", "2"), 
                   labels = c("Homdel", "Hetlos", "Diploid", "Gain", "Amp"), 
                   drop   = FALSE) +
  theme(legend.position   = "none")

# compare the variance among different groups using a nonparametri test
fligner.test(SCD~PTEN.CNV, data=PTEN.CNV.DNFA.RNAseq) # p-value = 0.0004647 (pten loss correlates with scd decrease?)
fligner.test(PTEN~PTEN.CNV, data=PTEN.CNV.DNFA.RNAseq) # p-value = p-value = 5.948e-08

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


