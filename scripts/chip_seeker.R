# TODO(dlroxe): replace 'source' line with 'library'
# line once everything is organized into a proper
# package.  Despite its location under 'scripts',
# this file assumes that its working directory
# is the project root directory (the parent directory
# of both 'scripts' and 'R')
source(file.path('R', 'ChIPseeker.R'))
chip_seeker_do_all()