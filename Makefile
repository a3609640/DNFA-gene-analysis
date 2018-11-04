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

$(starBinDir)/STAR: $(starSrcDir)
	cd $< && make 

$(starSrcDir): ${DNFA_starRoot}/2.6.0a.tar.gz
	cd $(<D) && tar -xzf $<

${DNFA_starRoot}/2.6.0a.tar.gz:
	$(MKDIR) $(@D)
	cd $(@D) && wget -nc https://github.com/alexdobin/STAR/archive/$(@F)

$(gtf): $$@.gz  # this requires .SECONDEXPANSION
	cd $(@D) && $(GUNZIP) $<

$(gtf).gz:
	$(MKDIR) $(@D)
	cd $(@D) && wget -nc $(ensemblBase)/gtf/homo_sapiens/$(@F)

$(fa): $$@.gz  # this requires .SECONDEXPANSION
	cd $(@D) && $(GUNZIP) $<

$(fa).gz:
	$(MKDIR) $(@D)
	cd $(@D) && wget -nc $(ensemblBase)/fasta/homo_sapiens/dna/$(@F)

# 'SA' is one of many alignment output files, but used as a proxy for the whole alignment step
reference_genome: ${DNFA_starRoot}/STARIndex/SA 

${DNFA_starRoot}/STARIndex/SA: $(starBinDir)/STAR $(fa) $(gtf) 
	$(MKDIR) $(@D)
	cd $(@D) && $(starBinDir)/STAR \
	  --runThreadN 12 \
	  --runMode genomeGenerate \
	  --genomeDir $(@D) \
	  --genomeFastaFiles $(word 2, $^) \
	  --sjdbGTFfile $(word 3, $^) \
	  --sjdbOverhang 75

clean: clean_star clean_refgenome

clean_star:
	rm -rf ${DNFA_starRoot}

clean_refgenome:
	rm -rf ${DNFA_generatedDataRoot}/referenceGenome

