install.packages('cgdsr')
library(cgdsr)

# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

test(mycgds)

# Get list of cancer studies at server
getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
allcancerstudies = getCancerStudies(mycgds)
# only retrieve tcga provisional data from server
tcga_provisional_studies <- allcancerstudies[ grep("(TCGA, Provisional)", allcancerstudies$name), ]

# get available case lists for a specific cancer study
# TODO please make a loop to list all cancer studies in the tcga_provisional_studies and their case lists
tcga_provisional_studies_caselist = getCaseLists(mycgds,tcga_provisional$cancer_study_id)[1,1]




skcm_tcga = getCancerStudies(mycgds)[194,1]
skcmcaselist = getCaseLists(mycgds,'skcm_tcga')
skcm_tcga_all = getCaseLists(mycgds,'skcm_tcga')[2,1]


# Get available genetic profiles
SKCMgeneticprofile = getGeneticProfiles(mycgds,'skcm_tcga')


# Get data slices for a specified list of genes, genetic profile and case list
BRAF.mutations=getProfileData(mycgds,c('BRAF'),
                              'skcm_tcga_mutations','skcm_tcga_all')
EGFR.mutations=getProfileData(mycgds,c('EGFR'),
                              'skcm_tcga_mutations','skcm_tcga_all')

SREBF1.CNV=getProfileData(mycgds,c('SREBF1'),
                          'skcm_tcga_linear_CNA','skcm_tcga_cna')
SREBF1.RNAseqZ=getProfileData(mycgds,c('ACACA','FASN','SCD','SREBF1'),
                              'skcm_tcga_rna_seq_v2_mrna_median_Zscores','skcm_tcga_all')


SREBF1.SCD.RNAseq=getProfileData(mycgds,c('ACACA','FASN','SCD','SREBF1','HMGCS1','HMGCR','LDLR','SREBF2','CD274','MITF'),
                                 'skcm_tcga_rna_seq_v2_mrna','skcm_tcga_all')
SREBF1.SCD.RPPA=getProfileData(mycgds,c('ACACA','FASN','SCD'),
                               'skcm_tcga_rppa','skcm_tcga_sequenced')
SREBF1.SCD.methy=getProfileData(mycgds,c('ACACA','FASN','SCD','SREBF1'),
                                'skcm_tcga_methylation_hm450','skcm_tcga_sequenced')
SREBF1.SCD.CNA=getProfileData(mycgds,c('ACACA','FASN','SCD','SREBF1'),
                              'skcm_tcga_linear_CNA','skcm_tcga_sequenced')

BRAF.RNAseq=cbind(BRAF.mutations, SREBF1.SCD.RNAseq)
EGFR.RNAseq=cbind(EGFR.mutations, SREBF1.SCD.RNAseq)

BRAF.RPPA=cbind(BRAF.mutations, SREBF1.SCD.RPPA)
BRAF.methy=cbind(BRAF.mutations, SREBF1.SCD.methy)
BRAF.CNA=cbind(BRAF.mutations, SREBF1.SCD.CNA)


levels(BRAF.RNAseq$BRAF) <- c(levels(BRAF.RNAseq$BRAF), "Mutated") 
levels(EGFR.RNAseq$EGFR) <- c(levels(EGFR.RNAseq$EGFR), "Mutated") 

levels(BRAF.RPPA$BRAF) <- c(levels(BRAF.RPPA$BRAF), "Mutated") 
levels(BRAF.methy$BRAF) <- c(levels(BRAF.methy$BRAF), "Mutated") 
levels(BRAF.CNA$BRAF) <- c(levels(BRAF.CNA$BRAF), "Mutated") 


BRAF.RNAseq$BRAF[BRAF.RNAseq$BRAF!='NA'] <- "Mutated"
EGFR.RNAseq$EGFR[EGFR.RNAseq$EGFR!='NA'] <- "Mutated"

BRAF.RPPA$BRAF[BRAF.RPPA$BRAF!='NA'] <- "Mutated"
BRAF.methy$BRAF[BRAF.methy$BRAF!='NA'] <- "Mutated"
BRAF.CNA$BRAF[BRAF.CNA$BRAF!='NA'] <- "Mutated"

p <- ggplot(BRAF.RNAseq, aes(x=BRAF, y=log2(SCD), color=BRAF))+ geom_violin(alpha = .5, trim=FALSE)+
    geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))

p <-ggplot(EGFR.RNAseq, aes(x=EGFR, y=log2(SCD), color=EGFR))+ geom_violin(alpha = .5, trim=FALSE)+
    geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))


cor.test(BRAF.RNAseq$CD274,BRAF.RNAseq$SCD, method = "pearson")


ggplot(BRAF.RPPA, aes(x=BRAF, y=SCD, fill=BRAF))+ geom_violin(alpha = .5, trim=FALSE)+
    geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))

ggplot(BRAF.methy, aes(x=BRAF, y=SREBF1, fill=BRAF))+ geom_violin(alpha = .5, trim=FALSE)+
    geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))

ggplot(BRAF.CNA, aes(x=BRAF, y=SCD, fill=BRAF))+ geom_boxplot()

SCD <- aggregate(SCD ~  BRAF, BRAF.SREBF1, median)
write.csv(BRAF.RPPA, "RPPA.csv")

SREBF1.mutation=getProfileData(mycgds,c('SREBF1','BRAF','NRAS','NF1'),
                               'skcm_tcga_mutations','skcm_tcga_sequenced')
levels(SREBF1.mutation$SREBF1) <- c(levels(SREBF1.mutation$SREBF1), "P783Q","L387M","G775W","S337Y","G89W","R1059S","Mutated") 
SREBF1.mutation["TCGA.D3.A2JG.06","SREBF1"]="P783Q"
SREBF1.mutation["TCGA.D3.A3C8.06","SREBF1"]="L387M"
SREBF1.mutation["TCGA.ER.A19B.06","SREBF1"]="G775W"
SREBF1.mutation["TCGA.ER.A19P.06","SREBF1"]="S337Y"
SREBF1.mutation["TCGA.FS.A1Z0.06","SREBF1"]="G89W"
SREBF1.mutation["TCGA.HR.A2OG.06","SREBF1"]="R1059S"

SREBF1.RNAseq=cbind(SREBF1.mutation, SREBF1.SCD.RNAseq)
colnames(SREBF1.RNAseq)[4] <- "SREBF1mut"
SREBF1.RNAseq=subset(SREBF1.RNAseq, !is.na(SREBF1mut))
SREBF1.RNAseq$SREBF1mut[SREBF1.RNAseq$SREBF1mut!='NaN'] <- "Mutated"


levels(SREBF1.RNAseq$BRAF) <- c(levels(SREBF1.RNAseq$BRAF), "Mutated") 
SREBF1.RNAseq$BRAF[SREBF1.RNAseq$BRAF!='NaN'] <- "Mutated"
SREBF1.BRAF.RNAseq=subset(SREBF1.RNAseq, SREBF1.RNAseq$BRAF=="Mutated")

ggplot(SREBF1.BRAF.RNAseq, aes(x=SREBF1mut, y=log2(FASN), color=SREBF1mut))+ geom_violin(alpha = .5, trim=FALSE)+
    geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))
SCD.medians <- aggregate(SCD ~  SREBF1mut, SREBF1.BRAF.RNAseq, median)

levels(SREBF1.RNAseq$NRAS) <- c(levels(SREBF1.RNAseq$NRAS), "Mutated")
SREBF1.RNAseq$NRAS[SREBF1.RNAseq$NRAS!='NaN'] <- "Mutated"

levels(SREBF1.RNAseq$NF1) <- c(levels(SREBF1.RNAseq$NF1), "Mutated")
SREBF1.RNAseq$NF1[SREBF1.RNAseq$NF1!='NaN'] <- "Mutated"






ggplot(SREBF1.RNAseq, aes(x=SREBF1mut, y=log2(SREBF1), color=SREBF1mut))+ geom_violin(alpha = .5, trim=FALSE)+
    geom_boxplot(alpha = .01, width = .1, position = position_dodge(width = .9))

write.csv(SREBF1.RNAseq, "SREBF1.RNAseq.csv")



# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,'skcm_tcga_sequenced')
SREBF1.clinc=cbind(SREBF1.mutation, myclinicaldata)
SREBF1.clinc=merge(SREBF1.mutation, myclinicaldata, all = TRUE)
write.csv(SREBF1.clinc, "SREBF1.clinc.csv")
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
