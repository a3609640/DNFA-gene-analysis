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
starSrcDir=${DNFA_starRoot}/STAR-2.6.1a/source
starBinDir=${DNFA_starRoot}/STAR-2.6.1a/bin

# FTP root for fetching human reference genome data from ENSEMBL
ensemblBase=ftp://ftp.ensembl.org/pub/release-94

# local output files where human gene annotation file (.gtf) and 
# human reference genome file (.fa)  from ENSEMBL will be stored
gtf=${DNFA_generatedDataRoot}/referenceGenome/Homo_sapiens.GRCh38.94.gtf
fa=${DNFA_generatedDataRoot}/referenceGenome/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# convenience definitions for commands with options (e.g. if permissions other than 755
# seem useful for 'mkdir' later on, it only needs to be changed in this one place)
GUNZIP=gunzip -k
MKDIR=mkdir -m 755 -p

# in accord with usual conventions, the first target is named 'all'
all: bamfiles

$(starBinDir)/STAR: | $(starSrcDir)
	cd $(word 1, $|) && $(MAKE) STAR STARlong install

$(starSrcDir): ${DNFA_starRoot}/2.6.1a.tar.gz
	cd $(<D) && tar -xzf $<

${DNFA_starRoot}/2.6.1a.tar.gz: | $$(@D)
	cd $(@D) && wget -nc https://github.com/alexdobin/STAR/archive/$(@F)

$(fa) $(gtf): $$@.gz
	cd $(@D) && $(GUNZIP) $<

$(gtf).gz: | $$(@D)
	cd $(@D) && wget -nc $(ensemblBase)/gtf/homo_sapiens/$(@F)

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
testIDs3 = 3_S4_L001 3_S4_L002 3_S4_L003 3_S4_L004
testIDs4 = 4_S6_L001 4_S6_L002 4_S6_L003 4_S6_L004
testIDs5 = 5_S5_L001 5_S5_L002 5_S5_L003 5_S5_L004
testIDs6 = 6_S1_L001 6_S1_L002 6_S1_L003 6_S1_L004
testIDs = $(testIDs1) $(testIDs2) $(testIDs3) $(testIDs4) $(testIDs5) $(testIDs6)

# TODO(dlroxe): probably, 'test1' and 'test2' should be respectively renamed
# to 'test1-40222594' and 'test2-40220751 .  That way they will match the
# unaltered contents of Testrun-pristine.tar.

# The 'bamfiles' target depends on a list of files with names
# like 'test1_S2_L003Aligned.sorted.bam'.
bamfiles: $(foreach id, $(testIDs), $(bamDir)/test$(id)Aligned.sorted.bam)
samfiles: \
 $(foreach id, $(testIDs1), $(samDir)/test1/test$(id)Aligned.out.sam) \
 $(foreach id, $(testIDs2), $(samDir)/test2/test$(id)Aligned.out.sam) \
 $(foreach id, $(testIDs3), $(samDir)/test3-40218815/test$(id)Aligned.out.sam) \
 $(foreach id, $(testIDs4), $(samDir)/test4-40218817/test$(id)Aligned.out.sam) \
 $(foreach id, $(testIDs5), $(samDir)/test5-40233238/test$(id)Aligned.out.sam) \
 $(foreach id, $(testIDs6), $(samDir)/test6-40220753/test$(id)Aligned.out.sam)

# For each .bam file we might want to build, locate a similarly-named
# .sam file and execute 'samtools' on it.

$(bamDir)/%.sorted.bam : $(samDir)/test1/%.out.sam | $$(@D)
	 $(sam2bam)

$(bamDir)/%.sorted.bam : $(samDir)/test2/%.out.sam | $$(@D)
	 $(sam2bam)

$(bamDir)/%.sorted.bam : $(samDir)/test3-40218815/%.out.sam | $$(@D)
	 $(sam2bam)

$(bamDir)/%.sorted.bam : $(samDir)/test4-40218817/%.out.sam | $$(@D)
	 $(sam2bam)

$(bamDir)/%.sorted.bam : $(samDir)/test5-40233238/%.out.sam | $$(@D)
	 $(sam2bam)

$(bamDir)/%.sorted.bam : $(samDir)/test6-40220753/%.out.sam | $$(@D)
	 $(sam2bam)

define sam2bam =
	 rm -f $@.tmp.* && samtools sort $< -o $@
endef

$(samDir)/%Aligned.out.sam : $(starBinDir)/STAR \
                             $(gtf) \
                             ${DNFA_raw_data_basedir}/%_R1_001.fastq.gz \
                             ${DNFA_raw_data_basedir}/%_R2_001.fastq.gz \
                             | reference_genome $$(@D)
	$< \
     --genomeDir ${DNFA_generatedDataRoot}/STARIndex \
     --sjdbGTFfile $(word 2, $^) \
     --runThreadN 12 \
     --outFileNamePrefix $(@D)/$(*F) \
     --readFilesIn $(word 3, $^) $(word 4, $^) \
     --readFilesCommand zcat \
     --twopassMode Basic

# Extract individual .fastq.gz files from parent .tar file, as needed.
${DNFA_raw_data_basedir}/%.fastq.gz:
	tar xvf ${DNFA_raw_data_tarfile} \
	  -C ${DNFA_raw_data_basedir} \
	  --transform="s/Testrun\/Nextseq\ Test-32743727\///;" \
	  Testrun/Nextseq\ Test-32743727/$(*D)/$(*F).fastq.gz

${DNFA_starRoot}:
	$(MKDIR) $@

${DNFA_generatedDataRoot}/referenceGenome ${DNFA_generatedDataRoot}/STARIndex:
	$(MKDIR) $@

$(samDir)/test1 $(samDir)/test2 $(samDir)/test3-40218815:
	$(MKDIR) $@

$(samDir)/test4-40218817 $(samDir)/test5-40233238 $(samDir)/test6-40220753:
	$(MKDIR) $@

$(bamDir):
	$(MKDIR) $@

# Use "make clean" to remove all generated output.
# Use sub-targets (e.g. "make clean_star") to remove selected
# generated output.
clean: clean_star clean_starindex clean_refgenome clean_samfiles clean_bamfiles

clean_star:
	rm -rf ${DNFA_starRoot}

clean_starindex:
	rm -rf ${DNFA_generatedDataRoot}/STARIndex

clean_refgenome:
	rm -rf ${DNFA_generatedDataRoot}/referenceGenome

clean_samfiles:
	rm -rf $(samDir)

clean_bamfiles:
	rm -rf $(bamDir)

.PHONY: all \
        clean clean_star clean_starindex clean_refgenome \
        clean_samfiles clean_bamfiles \
        reference_genome samfiles bamfiles
