# This Makefile is designed for use with GNU Make.
#
# See also https://www.gnu.org/software/make/manual/make.html .
#
# For a discussion of "automatic variables" such as $@, $<,
# $(@D), $(<F), and the like see in particular that section of the manual:
#
# https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html#Automatic-Variables
#

include configuration  # automatically grab environment variables for local env

.SECONDEXPANSION:  # to allow use of $@ in prequisite lists (escaped as $$@)

# output directories where STAR sources and binaries will be genereated
starSrcDir=${DNFA_starRoot}/STAR-2.6.1a/source
starBinDir=${DNFA_starRoot}/STAR-2.6.1a/bin

# FTP root for fetching human reference genome data from ENSEMBL
ensemblBase=ftp://ftp.ensembl.org/pub/release-94

# local output files where human gene annotation file (.gtf) and
# human reference genome file (.fa)  from ENSEMBL will be stored
gtf=${DNFA_generatedDataRoot}/referenceGenome/Homo_sapiens.GRCh38.94.gtf
fa=${DNFA_generatedDataRoot}/referenceGenome/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Stuff we might put under inst/extdata, but that is in general too
# big for RStudio to handle efficiently.  We'll make symlinks to this
# location instead.
extDataDir=${DNFA_generatedDataRoot}/r-extdata

# R project external data directory
projectExtDataDir=inst/extdata

# convenience definitions for commands with options (e.g. if permissions other than 755
# seem useful for 'mkdir' later on, it only needs to be changed in this one place)
GUNZIP=gunzip -k
MKDIR=mkdir -m 755 -p

# in accord with usual conventions, the first target is named 'all'
all: bamfiles $(extDataDir)/gene-counts.csv projectreadfiles

#############################################
# install the STAR aligner software package #
#############################################
# 'make ${DNFA_starRoot}/STAR-2.6.1a/bin/STAR' needs the source dir to exist
# (the '|' means it will not check the time at which the dir was modified)
# and will build STAR from the source dir fir STAR does not already exist.
$(starBinDir)/STAR: | $(starSrcDir)
	cd $(word 1, $|) && $(MAKE) STAR STARlong install

# To create the STAR source dir, we need the STAR .tar.gz file.  Once
# we have that, we go to ${DNFA_starRoot} and unpack the file.
$(starSrcDir): ${DNFA_starRoot}/2.6.1a.tar.gz
	cd $(<D) && tar -xzf $<

# To get the STAR .tar.gz file, first we need ${DNFA_starRoot} to
# exist (though again with '|', we do not care how recently it was
# modified).  Given the directory, we 'cd' there and 'wget' the
# file from github.
${DNFA_starRoot}/2.6.1a.tar.gz: | $$(@D)
	cd $(@D) && wget -nc https://github.com/alexdobin/STAR/archive/$(@F)

########################################################
# Prepare the reference genome file for STAR alignment #
########################################################
# If we want the Ensembl .fa or .gtf files, first we need their .gz
# forms.  Given those, we can go to the directory where we want the
# genome data to be, and unzip them.
$(fa) $(gtf): $$@.gz
	cd $(@D) && $(GUNZIP) $<

# To get the Ensembl .gtf.gz file, first we need the parent
# directory to exist (not caring when it last changed).
# Given that, we can 'cd' there and 'wget' the file.
$(gtf).gz: | $$(@D)
	cd $(@D) && wget -nc $(ensemblBase)/gtf/homo_sapiens/$(@F)

# To get the Ensembl .fa.gz file, first we need the parent
# directory to exist (not caring when it last changed).
# Given that, we can 'cd' there and 'wget' the file.
$(fa).gz: | $$(@D)
	cd $(@D) && wget -nc $(ensemblBase)/fasta/homo_sapiens/dna/$(@F)

#######################
# generate STAR Index #
#######################
# 'SA' (containing a suffix array) is one of many alignment output files,
# but used as a proxy for the whole alignment step
reference_genome: ${DNFA_generatedDataRoot}/STARIndex/SA

# If we want the SA file, first STAR and the Ensembl .gtf and .fa files must
# exist.  If we need any of those, we make them using the rules above.  If
# any of them changed more recently than the SA file, we assume there's
# something new and rebuild the SA file.  We also need the STARIndex dir to
# exist, but we don't care how recently it changed.  Given all the prereqs,
# we can 'cd' to STARIndex and execute STAR to build the reference genome.

# Notes about the parameters
# --sjdbOverhang species the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads. We used 2x76b paired-end reads, the ideal value is 76-1=75.

${DNFA_generatedDataRoot}/STARIndex/SA: $(starBinDir)/STAR $(fa) $(gtf) | $$(@D)
	cd $(@D) && $< \
	  --runThreadN 12 \
	  --runMode genomeGenerate \
	  --genomeDir $(@D) \
	  --genomeFastaFiles $(word 2, $^) \
	  --sjdbGTFfile $(word 3, $^) \
	  --sjdbOverhang 75

# The next few lines define variables related to .sam and .bam files.
# The rules for building these files will be clever with pattern
# matching, but one important thing will be for the Makefile to know
# exactly what list of files it's trying to build.  These variables
# define output locations and unique test identifiers.
testIDs1 = 1_S2_L001 1_S2_L002 1_S2_L003 1_S2_L004
testIDs2 = 2_S3_L001 2_S3_L002 2_S3_L003 2_S3_L004
testIDs3 = 3_S4_L001 3_S4_L002 3_S4_L003 3_S4_L004
testIDs4 = 4_S6_L001 4_S6_L002 4_S6_L003 4_S6_L004
testIDs5 = 5_S5_L001 5_S5_L002 5_S5_L003 5_S5_L004
testIDs6 = 6_S1_L001 6_S1_L002 6_S1_L003 6_S1_L004
testIDs = $(testIDs1) $(testIDs2) $(testIDs3) $(testIDs4) $(testIDs5) $(testIDs6)

# The 'bamfiles' target depends on a list of files with names
# like 'test1_S2_L003Aligned.sorted.bam'.

# The full lists of .sam and .bam files are generated here with 'foreach'
# loops that match up test IDs (from above) with the test directories that
# contain those IDs in the input Testrun-pristine.tar file.  If we want all
# the .sam files, we want each individual .sam file generated by all the
# 'foreach' loops below.

bamfiles: $(foreach id, $(testIDs), $(extDataDir)/test$(id)Aligned.sorted.bam)
samfiles: $(foreach id, $(testIDs), $(extDataDir)/test$(id)Aligned.out.sam)
readfiles: $(foreach id, $(testIDs), $(extDataDir)/test$(id)ReadsPerGene.out.tab)
projectreadfiles: $(foreach id, $(testIDs), $(projectExtDataDir)/test$(id)ReadsPerGene.out.tab)

# Likewise, the full list of .fastq.gz files is generated here.  This also
# permits all the .fastq.gz files to be extracted from the archive file with
# 'make fastqfiles', and furthermore makes them explicit targets.  If these were
# only intermediate targets, 'make' might delete them after completing work
# on the .sam files that need them.
fastqfiles: \
 $(foreach n, 1 2, \
  $(foreach id, $(testIDs1), ${DNFA_raw_data_basedir}/test1-40222594/test$(id)_R$(n)_001.fastq.gz) \
  $(foreach id, $(testIDs2), ${DNFA_raw_data_basedir}/test2-40220751/test$(id)_R$(n)_001.fastq.gz) \
  $(foreach id, $(testIDs3), ${DNFA_raw_data_basedir}/test3-40218815/test$(id)_R$(n)_001.fastq.gz) \
  $(foreach id, $(testIDs4), ${DNFA_raw_data_basedir}/test4-40218817/test$(id)_R$(n)_001.fastq.gz) \
  $(foreach id, $(testIDs5), ${DNFA_raw_data_basedir}/test5-40233238/test$(id)_R$(n)_001.fastq.gz) \
  $(foreach id, $(testIDs6), ${DNFA_raw_data_basedir}/test6-40220753/test$(id)_R$(n)_001.fastq.gz) )

#################################################
# sort SAM files from the STAR alignment output #
#################################################
# For each .bam file we might want to build, locate a similarly-named
# .sam file and execute 'samtools' on it.  A given .bam file we want to
# build could require a .sam file from any of the test directories.
#
# So, we repeat the sam->bam target rule for each such directory, but
# we let each rule share the same $(sam2bam) recipe.
#
# Note the pattern matching.  For example, "%" could match "test2_S3_L001".
# In that case, where "%" has to mean the same string for both $.sorted.bam
# and $.out.sam, the second rule below (with 'test2' in its name) would
# match and cause the recipe to be invoked.

$(extDataDir)/%.sorted.bam : $(extDataDir)/%.out.sam | $$(@D)
	 rm -f $@.tmp.* && samtools sort $< -o $@

# To make a particular .sam file, first we need STAR, the Ensembl .gtf file,
# and two corresponding .fastq.gz files.  If any of those were updated more
# recently than our .sam file, we will remake the .sam file.  We also need
# a parent directory and a reference genome, but we don't care how recently
# they were changed.  Given all those prereqs, we can invoke STAR to read
# the inputs and generate the output in the location we choose.
#
# Again note the pattern matching here.  In this case, the wildcard '*'
# will end up matching a string like "test3-40218815", and the pattern
# matcher "%" will end up matching a string (for both .sam and .fastq.gz)
# like "test3_S4_L002"

#############################################
# Perform STAR alignment on the FASTQ files #
#############################################
# parameters used for STAR alignment #
# --readFilesCommand zcat option inputs readfiles that are compressed as gzipped files (*.gz))
# default output files include Aligned.out.sam (alignments in standard SAM format)
# --quantMode TranscriptomeSAM optionwill output alignments translated into transcript coordinates in the Aligned.toTranscriptome.out.bamfile (in addition to alignments in genomic coordinates in Aligned.*.sam
# --quantMode GeneCounts option will count number reads per gene while mapping.
# --quantMode TranscriptomeSAM GeneCounts option will get both the Aligned.toTranscriptome.out.bam and ReadsPerGene.out.tab output files.

$(extDataDir)/%ReadsPerGene.out.tab $(extDataDir)/%Aligned.out.sam : $(starBinDir)/STAR \
                             $(gtf) \
                             ${DNFA_raw_data_basedir}/*/%_R1_001.fastq.gz \
                             ${DNFA_raw_data_basedir}/*/%_R2_001.fastq.gz \
                             | reference_genome $$(@D)
	$< \
     --genomeDir ${DNFA_generatedDataRoot}/STARIndex \
     --sjdbGTFfile $(word 2, $^) \
     --runThreadN 12 \
     --outFileNamePrefix $(@D)/$(*F) \
     --quantMode TranscriptomeSAM GeneCounts \
     --readFilesIn $(word 3, $^) $(word 4, $^) \
     --readFilesCommand zcat \
     --twopassMode Basic

$(projectExtDataDir)/%ReadsPerGene.out.tab : $(extDataDir)/%ReadsPerGene.out.tab
	ln -s $< $@

$(extDataDir)/gene-counts.csv: R/DESeqDataPreparation.R readfiles
	Rscript $<

# To make a given .fastq.gz file, first we need the .tar file that contains
# it, and a directory in which to put it.  Given that, we can execute a tar
# command to perform the extraction.  Note again the pattern matching, which
# in this case will turn out to match a string like
# "test6-40220753/test6_S1_L004_R1_001".  (The '%' matcher can in general
# match patterns that span directories.)
#
# Mark the extracted file (with 'touch') as more recently modified than the
# archive, prevent it from being re-extracted (with an old 'modified' date)
# in every run.
${DNFA_raw_data_basedir}/%.fastq.gz: ${DNFA_raw_data_tarfile} | ${DNFA_raw_data_basedir}
	tar -xvf $< \
	  -C ${DNFA_raw_data_basedir} \
	  --transform="s/Testrun\/Nextseq\ Test-32743727\///;" \
	  Testrun/Nextseq\ Test-32743727/$(*D)/$(*F).fastq.gz \
	&& touch $@

# For each of the following targets, if we want it to exist we create
# it as a directory with 'mkdir'.
${DNFA_starRoot}:
	$(MKDIR) $@

# The trailing '/' is super-important here to force pattern matching
# to stop before reaching a filename.  Otherwise, 'make' will create
# a *directory* named ....fastq.gz, instead of creating a parent dir
# into which that file can be extracted.
${DNFA_raw_data_basedir}/test%/: | ${DNFA_raw_data_basedir}
	$(MKDIR) $@

${DNFA_raw_data_basedir}:
	$(MKDIR) $@

${DNFA_generatedDataRoot}/referenceGenome ${DNFA_generatedDataRoot}/STARIndex:
	$(MKDIR) $@

$(extDataDir):
	$(MKDIR) $@

# Use "make clean" to remove all generated output.
# Use sub-targets (e.g. "make clean_star") to remove selected
# generated output.
clean: clean_starindex clean_refgenome clean_samfiles clean_bamfiles

# Oops -- that's a dangerous directive that will delete the whole system
# if run with root permissions and the environment variable is not set.
#
#
#clean_star:
#	rm -rf ${DNFA_starRoot}

clean_starindex:
	rm -rf ${DNFA_generatedDataRoot}/STARIndex

clean_refgenome:
	rm -rf ${DNFA_generatedDataRoot}/referenceGenome

clean_samfiles:
	rm $(extDataDir)/*.sam

clean_bamfiles:
	rm $(extDataDir)/*.bam

clean_fastqfiles:
	rm -rf ${DNFA_raw_data_basedir}

# The .PHONY target is special -- see also the GNU Makefile documentation.
# Basically it helps to disambiguate the case where we have a target named
# 'clean' but there may be a file named 'clean' that could be interpreted
# as a target.  This makes it clear that the 'clean' target is not meant
# to correspond to a filename.
.PHONY: all \
        clean clean_starindex clean_refgenome \
        clean_samfiles clean_bamfiles \
        reference_genome samfiles bamfiles
