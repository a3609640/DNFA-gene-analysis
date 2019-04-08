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

##################################################################
## plot DNFA gene expression across all TCGA provisional groups ##
##################################################################
plot.DNFA.provisional.tcga <- function(DNFA) {
  ### retrieve all TCGA provisional study groups
  tcga.pro.studies <- getCancerStudies(mycgds)[
    grep("(TCGA, Provisional)", getCancerStudies(mycgds)$name), ]
  ### "tcag_study_list" contains all the tcga cancer studies
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
      grep("tcga_rna_seq_v2_mrna",
           tcga.pro.geneticprofile[[x]]$genetic_profile_id), ][1, 1]
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
  df2 <- na.omit(df2)
  # plot DNFA gene expression across all TCGA groups ##
  m <- paste0(DNFA, ".", DNFA)
  mean <- within(df2[df2$DNFAgene == m,], # TCGAstudy is one column in df2
                 TCGAstudy <- reorder(TCGAstudy, log2(RNAseq), median))
  a <- levels(mean$TCGAstudy)
  data.long <- mean
  data.long$names <- rownames(data.long)
  data.wide <- dcast(data.long,  
                     DNFAgene + names ~ TCGAstudy, 
                     value.var = "RNAseq")
  ### write.csv(data.wide, file = paste(DNFA, ".csv", sep = ""))
  ### with highlight on skin cancer
  colors <- ifelse(a == "skcm_tcga", "red", "black")
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
           y = paste0("log2(", DNFA, " RNA counts)")) +
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

plot.DNFA.provisional.tcga("SCD")
lapply(DNFA.gene, plot.DNFA.provisional.tcga)

################################################################
## plot DNFA gene expression across all TCGA pancancer groups ##
################################################################
plot.DNFA.pan.tcga <- function(DNFA){
  # Get DNFA RNAseq data from all TCGA study groups
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
  df2 <- na.omit(df2)
  ### plot DNFA gene expression across all TCGA groups ##
  m <- paste0(DNFA, ".", DNFA)
  mean <- within(df2[df2$DNFAgene == m,], # TCGAstudy is one column in df2
                 TCGAstudy <- reorder(TCGAstudy, log2(RNAseq), median))
  a <- levels(mean$TCGAstudy)
  ### write.csv(mean, file = paste(DNFA, ".mean", sep=""))
  ### with highlight on skin cancer
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
           y = paste0("log2(", DNFA, " RNA counts)")) +
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
  bug.data.set <- names(pro.tcga.study.list) %in% c("meso_tcga",
    "pcpg_tcga",
    "ucs_tcga")
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
  number <- round(nrow(df)/5)
  df$Group[df$RNAseq < quantile(df$RNAseq, prob = 0.2)] = "Bottom 20%"
  df$Group[df$RNAseq > quantile(df$RNAseq, prob = 0.8)] = "Top 20%"
  df$SurvObj <- with(df, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))
  #  df <- na.omit(df)
  km <- survfit(SurvObj ~ df$Group, data = df, conf.type = "log-log")
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 0)
  p.val <- 1 - pchisq(stats$chisq, length(stats$n) - 1)
  p.val <- signif(p.val, 3)
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
        legend.justification = c(1,1),
        legend.position      = c(1,1))+
      guides(fill = FALSE) +
      scale_color_manual(values = c("red", "blue"),
        name   = paste(DNFA, "mRNA expression"),
        breaks = c("Bottom 20%", "Top 20%"),
        labels = c(paste("Bottom 20%, n =", number),
          paste("Top 20%, n =", number))) +
      geom_point(size = 0.25) +
      annotate("text",
        x     = 250,
        y     = 0.80,
        label = paste("log-rank test, p.val = ", p.val),
        size  = 4.5,
        hjust = 1,
        fontface = "bold"))
  
  # rho = 1 the Gehan-Wilcoxon test
  stats <- survdiff(SurvObj ~ df$Group, data = df, rho = 1)
  print(DNFA)
  print(stats)
}

plot.km.all.tcga("SCD")
sapply(DNFA.gene, plot.km.all.tcga)





