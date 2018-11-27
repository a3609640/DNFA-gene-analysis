install.packages('cgdsr')
library(cgdsr)
library(ggplot2)
library(reshape2)

# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
test(mycgds)
# Get list of cancer studies at server
getCancerStudies(mycgds)
# Get cases from TCGA provisional studies only
tcga_provisional_studies = getCancerStudies(mycgds)[grep("(TCGA, Provisional)", 
                                                         getCancerStudies(mycgds)$name), ]
# "tcag_study_list" is a vector containing all the tcga cancer studies that I would to analyze for DNFA gene expression
tcag_study_list = tcga_provisional_studies$cancer_study_id
names(tcag_study_list) <- tcag_study_list

caselist = function (x) getCaseLists(mycgds, x)
geneticprofile = function (x) getGeneticProfiles(mycgds, x)
# use lappy to pull out all the caselists within tcag_study_list
# lappy will return a large list, each element in that list is a data-table
tcag_provisional_caselist = lapply (tcag_study_list, caselist)
tcag_provisional_geneticprofile = lapply (tcag_study_list, geneticprofile)
# for example, tcag_provisional_caselist[[1]] shows the dataframe of caselist in laml study group.
# we want to choose case_list_id that is labeled with laml_tcga_rna_seq_v2_mrna, tcag_provisional_caselist[[1][8,1]
  # a =tcag_provisional_caselist[[1]][grep("tcga_rna_seq_v2_mrna",
  #                                       tcag_provisional_caselist[[1]]$case_list_id), ][1,1]
  # b =tcag_provisional_geneticprofile[[1]][grep("mRNA expression \\(RNA Seq V2 RSEM\\)", 
  # double backslash \\ suppress the special meaning of ( ) in regular expression
  #                                             tcag_provisional_geneticprofile[[1]]$genetic_profile_name), ][1,1]
# how do we do this for all study groups from [[1]] to  [[32]]?
caselist_RNAseq = function(x)
{tcag_provisional_caselist[[x]][grep("tcga_rna_seq_v2_mrna",
                                       tcag_provisional_caselist[[x]]$case_list_id), ][1,1]}
geneticprofile_RNAseq = function(x)
  {tcag_provisional_geneticprofile[[x]][grep("mRNA expression \\(RNA Seq V2 RSEM\\)",
                                             tcag_provisional_geneticprofile[[x]]$genetic_profile_name), ][1,1]}
# test the functions: caselist_RNAseq () and geneticprofile_RNAseq ()
    # caselist_RNAseq = caselist_RNAseq ('acc_tcga')
    # geneticprofile_RNAseq = geneticprofile_RNAseq ('acc_tcga')
# We wrap two functions: geneticprofile_RNAseq(x), caselist_RNAseq(x) within TCGA_ProfileData_RNAseq(x)
TCGA_ProfileData_RNAseq = function(geneticprofile, caselist) {getProfileData(mycgds, 'FASN', geneticprofile, caselist)}
TCGA_RNAseq = function(x) {TCGA_ProfileData_RNAseq(geneticprofile_RNAseq(x), caselist_RNAseq(x))}  
# test the wrapping function
  # data5 = TCGA_RNAseq ('acc_tcga')
TCGA_RNAseq_all_Studies = lapply (tcag_study_list, TCGA_RNAseq)

# use the melt function from reshape2 package. 
df <- melt(TCGA_RNAseq_all_Studies)
# Separate boxplots for each data.frame
qplot(factor(L1), log2(value), data = df, geom = "boxplot")
# violin graph for FASN gene expression across all tumor types, and order the X axis based on the gene expression level
mean <- within(df, L1 <-  reorder(L1, log2(value), median))
ggplot(mean, aes(x=L1, y=log2(value), fill=L1)) + 
  geom_violin(alpha = .5, trim=FALSE)+
  geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))+
  theme_bw()+
  labs(x = "Tumor types (TCGA)",
       y = paste("log2(FASN RNA counts)")) +
  theme(axis.title=element_text(face="bold",size=9,color="black"),
        axis.text=element_text(size=9,angle = 45, hjust = 1, face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 9, colour = "black"),
        legend.position = "none")

#################################################################################################
## get DNFA gene expression from SKCM group and check the correlation with oncogenic mutations ##
#################################################################################################
skcm_case = getCaseLists(mycgds,'skcm_tcga')
skcm_case = getCaseLists(mycgds,'skcm_tcga')
skcm_tcga_all = getCaseLists(mycgds,'skcm_tcga')[2,1]

# Get available genetic profiles
SKCMgeneticprofile = getGeneticProfiles(mycgds,'skcm_tcga')
DNFA.RNAseq=getProfileData(mycgds,c('ACACA','FASN','SCD','ACLY','ACSS2','ACSL1','LDLR','SREBF1','SREBF2','MITF'),
                           'skcm_tcga_rna_seq_v2_mrna',
                           'skcm_tcga_all')

# Get data slices for a specified list of genes, genetic profile and case list
BRAF.mutations=getProfileData(mycgds, 'BRAF',
                              'skcm_tcga_mutations',
                              'skcm_tcga_all')

NRAS.mutations=getProfileData(mycgds, 'NRAS',
                              'skcm_tcga_mutations',
                              'skcm_tcga_all')

AKT.mutations=getProfileData(mycgds, 'AKT',
                              'skcm_tcga_mutations',
                              'skcm_tcga_all')

TP53.mutations=getProfileData(mycgds, 'TP53',
                             'skcm_tcga_mutations',
                             'skcm_tcga_all')

BRAF.CNV=getProfileData(mycgds, 'BRAF',
                        "skcm_tcga_gistic",
                        "skcm_tcga_all")

NRAS.CNV=getProfileData(mycgds, "NRAS",
                        "skcm_tcga_gistic",
                        "skcm_tcga_all")

PTEN.CNV=getProfileData(mycgds, "PTEN",
                        "skcm_tcga_gistic",
                        "skcm_tcga_all")

BRAF.mutations.DNFA.RNAseq=cbind(BRAF.mutations, DNFA.RNAseq)
NRAS.mutations.DNFA.RNAseq=cbind(NRAS.mutations, DNFA.RNAseq)
AKT.mutations.DNFA.RNAseq=cbind(AKT.mutations, DNFA.RNAseq)
TP53.mutations.DNFA.RNAseq=cbind(TP53.mutations, DNFA.RNAseq)

BRAF.CNV.DNFA.RNAseq=cbind(BRAF.CNV, DNFA.RNAseq)
NRAS.CNV.DNFA.RNAseq=cbind(NRAS.CNV, DNFA.RNAseq)
PTEN.CNV.DNFA.RNAseq=cbind(PTEN.CNV, DNFA.RNAseq)

####################################################################################################################
####################################################################################################################
levels(BRAF.mutations.DNFA.RNAseq$BRAF) <- c(levels(BRAF.mutations.DNFA.RNAseq$BRAF), "Mutated", "Wildtype")
BRAF.mutations.DNFA.RNAseq$BRAF[BRAF.mutations.DNFA.RNAseq$BRAF!='NA'] <- "Mutated"
BRAF.mutations.DNFA.RNAseq$BRAF[is.na(BRAF.mutations.DNFA.RNAseq$BRAF)] <- "Wildtype" 
levels(droplevels(BRAF.mutations.DNFA.RNAseq$BRAF))
ggplot(BRAF.mutations.DNFA.RNAseq, aes(x=BRAF, y=log2(SCD), color=BRAF)) + 
  geom_violin(alpha = .5, trim=FALSE) + 
  geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))

levels(NRAS.mutations.DNFA.RNAseq$NRAS) <- c(levels(NRAS.mutations.DNFA.RNAseq$NRAS), "Mutated", "Wildtype")
NRAS.mutations.DNFA.RNAseq$NRAS[NRAS.mutations.DNFA.RNAseq$NRAS!='NA'] <- "Mutated"
NRAS.mutations.DNFA.RNAseq$NRAS[is.na(NRAS.mutations.DNFA.RNAseq$NRAS)] <- "Wildtype" 
ggplot(NRAS.mutations.DNFA.RNAseq, aes(x=NRAS, y=log2(SCD), color=NRAS)) + 
  geom_violin(alpha = .5, trim=FALSE) + 
  geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))

levels(AKT.mutations.DNFA.RNAseq$AKT) <- c(levels(AKT.mutations.DNFA.RNAseq$AKT), "Mutated", "Wildtype")
AKT.mutations.DNFA.RNAseq$AKT[AKT.mutations.DNFA.RNAseq$AKT!='NA'] <- "Mutated"
AKT.mutations.DNFA.RNAseq$AKT[is.na(AKT.mutations.DNFA.RNAseq$AKT)] <- "Wildtype" 
ggplot(AKT.mutations.DNFA.RNAseq, aes(x=AKT, y=log2(SCD), color=AKT)) + 
  geom_violin(alpha = .5, trim=FALSE)+
  geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))

levels(TP53.mutations.DNFA.RNAseq$TP53) <- c(levels(TP53.mutations.DNFA.RNAseq$TP53), "Mutated", "Wildtype")
TP53.mutations.DNFA.RNAseq$TP53[TP53.mutations.DNFA.RNAseq$TP53!='NA'] <- "Mutated"
TP53.mutations.DNFA.RNAseq$TP53[is.na(TP53.mutations.DNFA.RNAseq$TP53)] <- "Wildtype" 
ggplot(TP53.mutations.DNFA.RNAseq, aes(x=TP53, y=log2(SCD), color=TP53)) + 
  geom_violin(alpha = .5, trim=FALSE) +
  geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))

####################################################################################################################
####################################################################################################################
BRAF.CNV.DNFA.RNAseq$BRAF=as.factor(BRAF.CNV.DNFA.RNAseq$BRAF)
ggplot(BRAF.CNV.DNFA.RNAseq, aes(x=BRAF, y=log2(SCD), color=BRAF)) + 
  geom_violin(alpha = .5, trim=FALSE)+
  geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))

NRAS.CNV.DNFA.RNAseq$NRAS=as.factor(NRAS.CNV.DNFA.RNAseq$NRAS)
ggplot(NRAS.CNV.DNFA.RNAseq, aes(x=NRAS, y=log2(SCD), color=NRAS)) + 
  geom_violin(alpha = .5, trim=FALSE) +
  geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))

PTEN.CNV.DNFA.RNAseq$PTEN=as.factor(PTEN.CNV.DNFA.RNAseq$PTEN)
ggplot(PTEN.CNV.DNFA.RNAseq, aes(x=PTEN, y=log2(SCD), color=PTEN)) + 
  geom_violin(alpha = .5, trim=FALSE) +
  geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))

####################################################################################################################
####################################################################################################################

# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,'skcm_tcga_sequenced')
SREBF1.clinc=cbind(SREBF1.mutation, myclinicaldata)
SREBF1.clinc=merge(SREBF1.mutation, myclinicaldata, all = TRUE)
write.csv(SREBF1.clinc, "SREBF1.clinc.csv")
```


