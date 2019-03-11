# TODO(dlroxe): replace 'source' line with 'library'
# line once everything is organized into a proper
# package.  Despite its location under 'scripts',
# this file assumes that its working directory
# is the project root directory (the parent directory
# of both 'scripts' and 'R')

rm(list=ls())  # guarantee from-scratch execution: clear the global environment
source(file.path('R', 'DESeqandFactoMineRforASOSpecificity.R'))

analyze_aso_specificity()
