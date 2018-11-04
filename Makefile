# This Makefile is designed for use with GNU Make.
#
# See also https://www.gnu.org/software/make/manual/make.html .
#
# For a discussion of "automatic variables" such as $@, $<,
# $(@D), $(<F), and the like see in particular that section of the manual:
#
# https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html#Automatic-Variables
#

.SECONDEXPANSION:  # to allow use of $@ in prequisite lists (escaped as $$@)

# output directories where STAR sources and binaries will be genereated
starSrcDir=$(DNFA_starRoot)/STAR-2.6.0a/source
starBinDir=$(DNFA_starRoot)/STAR-2.6.0a/bin/Linux_x86_64

# FTP root for fetching human genome data from ENSEMBL
ensemblBase=ftp://ftp.ensembl.org/pub/release-94

# local output files where human genome GTF and FA files from ENSEMBL will be stored
gtf=${DNFA_generatedDataRoot}/referenceGenome/Homo_sapiens.GRCh38.94.gtf
fa=${DNFA_generatedDataRoot}/referenceGenome/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# convenience definitions for commands with options (e.g. if permissions other than 755
# seem useful for 'mkdir' later on, it only needs to be changed in this one place)
GUNZIP=gunzip -k
MKDIR=mkdir -m 755 -p

# in accord with usual conventions, the first target is named 'all'
all: reference_genome

$(starBinDir)/STAR: | $(starSrcDir)
	cd $(word 1, $|) && $(MAKE)

$(starSrcDir): ${DNFA_starRoot}/2.6.0a.tar.gz
	cd $(<D) && tar -xzf $<

${DNFA_starRoot}/2.6.0a.tar.gz: | $$(@D)
	cd $(@D) && wget -nc https://github.com/alexdobin/STAR/archive/$(@F)

$(gtf): $$@.gz
	cd $(@D) && $(GUNZIP) $<

$(gtf).gz: | $$(@D)
	cd $(@D) && wget -nc $(ensemblBase)/gtf/homo_sapiens/$(@F)

$(fa): $$@.gz
	cd $(@D) && $(GUNZIP) $<

$(fa).gz: | $$(@D)
	cd $(@D) && wget -nc $(ensemblBase)/fasta/homo_sapiens/dna/$(@F)

# 'SA' (containing a suffix array) is one of many alignment output files,
# but used as a proxy for the whole alignment step
reference_genome: ${DNFA_generatedDataRoot}/STARIndex/SA 

${DNFA_generatedDataRoot}/STARIndex/SA: $(starBinDir)/STAR $(fa) $(gtf) | $$(@D)
	cd $(@D) && $< \
	  --runThreadN 12 \
	  --runMode genomeGenerate \
	  --genomeDir $(@D) \
	  --genomeFastaFiles $(word 2, $^) \
	  --sjdbGTFfile $(word 3, $^) \
	  --sjdbOverhang 75

samDir=${DNFA_generatedDataRoot}/Analysis/2-pass
bamDir=${DNFA_generatedDataRoot}/Analysis/Samsort
testIDs1 = 1_S2_L001 1_S2_L002 1_S2_L003 1_S2_L004
testIDs2 = 2_S3_L001 2_S3_L002 2_S3_L003 2_S3_L004
testIDs = $(testIDs1) $(testIDs2)

# The 'bamfiles' target depends on a list of files with names
# like 'test1_S2_L003Aligned.sorted.bam'.
bamfiles: $(foreach id, $(testIDs), $(bamDir)/test$(id)Aligned.sorted.bam)
samfiles: samfiles1 samfiles2
samfiles1: $(foreach id, $(testIDs1), $(samDir)/test1/test$(id)Aligned.out.sam)
samfiles2: $(foreach id, $(testIDs2), $(samDir)/test2/test$(id)Aligned.out.sam)

# For each .bam file we might want to build, locate a similarly-named
# .sam file and execute 'samtools' on it.

$(bamDir)/%.sorted.bam : $(samDir)/test1/%.out.sam | $$(@D)
	 $(sam2bam)

$(bamDir)/%.sorted.bam : $(samDir)/test2/%.out.sam | $$(@D)
	 $(sam2bam)

define sam2bam:
	 rm -f $@.tmp.* && samtools sort $< -o $@
endef

$(samDir)/test1/%Aligned.out.sam : $(starBinDir)/STAR \
                                   ${DNFA_raw_data_basedir}/test1/%_R1_001.fastq.gz \
                                   ${DNFA_raw_data_basedir}/test1/%_R2_001.fastq.gz \
                                   | $(gtf) reference_genome $$(@D)
	$(alignTestData)

$(samDir)/test2/%Aligned.out.sam : $(starBinDir)/STAR \
                                   ${DNFA_raw_data_basedir}/test2/%_R1_001.fastq.gz \
                                   ${DNFA_raw_data_basedir}/test2/%_R2_001.fastq.gz \
                                   | $(gtf) reference_genome $$(@D)
	$(alignTestData)

define alignTestData =
	$< \
     --genomeDir ${DNFA_generatedDataRoot}/STARIndex \
     --sjdbGTFfile $(word 1, $|) \
     --runThreadN 12 \
     --outFileNamePrefix $(@D)/$(*F) \
     --readFilesIn $(word 2, $^) $(word 3, $^) \
     --readFilesCommand zcat \
     --twopassMode Basic
endef

${DNFA_starRoot}:
	$(MKDIR) $@

${DNFA_generatedDataRoot}/referenceGenome ${DNFA_generatedDataRoot}/STARIndex:
	$(MKDIR) $@

$(samDir)/test1 $(samDir)/test2 $(bamDir):
	$(MKDIR) $@

# Use "make clean" to remove all generated output.
# Use sub-targets (e.g. "make clean_star") to remove selected
# generated output.
clean: clean_star clean_starindex clean_refgenome clean_bamfiles

clean_star:
	rm -rf ${DNFA_starRoot}

clean_starindex:
	rm -rf ${DNFA_generatedDataRoot}/STARIndex

clean_refgenome:
	rm -rf ${DNFA_generatedDataRoot}/referenceGenome

clean_bamfiles:
	rm -rf $(bamDir)

.PHONY: all clean clean_star clean_starindex clean_refgenome clean_bamfiles samfiles bamfiles
