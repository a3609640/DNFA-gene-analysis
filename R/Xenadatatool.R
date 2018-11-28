# Download UCSC Xena Datasets and load them into R by UCSCXenaTools
library(UCSCXenaTools)
data(XenaData)
# XenaGenerate function, which generate XenaHub object from XenaData data frame.
head(XenaData)
XenaGenerate()
xe <- XenaGenerate()
# explore hosts(), cohorts() and datasets() inside XenaHub object
hosts(xe)
head(cohorts(xe))
head(datasets(xe))

# getTCGAdata provide an easy way to download TCGA datasets
args(getTCGAdata)
# See all available project id, please use
availTCGA("ProjectID")
# mRNASeq: logical. if TRUE, download mRNASeq data
# mRNASeqType: character vector.
# can be one, two or three in c("normalized", "pancan normalized", "percentile")
# download files to directory /tmp/RtmpyOXWWh.
# /tmp/RtmpyOXWWh/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz,
PANACAN <- getTCGAdata(project = "PANCAN", clinical = FALSE,
            mRNASeq = TRUE, mRNAArray = FALSE,
            mRNASeqType = "normalized", miRNASeq = FALSE,
            exonRNASeq = FALSE, RPPAArray = FALSE,
            Methylation = FALSE, MethylationType = "450K",
            GeneMutation = FALSE, SomaticMutation = FALSE,
            download = TRUE)
