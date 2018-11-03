starSrcDir=$(DNFA_starRoot)/STAR-2.6.0a/source
starBinDir=$(DNFA_starRoot)/STAR-2.6.0a/bin/Linux_x86_64

ensemblBase=ftp://ftp.ensembl.org/pub/release-94
genomeBase=$(ensemblBase)/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz
genomeAssembly=$(ensemblBase)/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

refGeneDir=${DNFA_generatedDataRoot}/referenceGenome

all: reference_genome

$(starBinDir)/STAR: $(starSrcDir)
	cd $< && make 

$(starSrcDir): ${DNFA_starRoot}/2.6.0a.tar.gz
	$(shell cd $(DNFA_starRoot) && tar -xzf $<)

${DNFA_starRoot}/2.6.0a.tar.gz:
	$(shell mkdir -m 755 -p ${DNFA_starRoot} && cd ${DNFA_starRoot} && wget -nc https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz)

$(refGeneDir)/Homo_sapiens.GRCh38.94.gtf : $(refGeneDir)/Homo_sapiens.GRCh38.94.gtf.gz
	$(shell cd $(refGeneDir) && gunzip -k $<)

$(refGeneDir)/Homo_sapiens.GRCh38.94.gtf.gz:
	$(shell mkdir -m 755 -p $(refGeneDir) && cd $(refGeneDir) && wget -nc ${genomeBase})

$(refGeneDir)/Homo_sapiens.GRCh38.dna.primary_assembly.fa : $(refGeneDir)/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	$(shell cd $(refGeneDir) && gunzip -k $<)

$(refGeneDir)/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz :
	$(shell mkdir -m 755 -p $(refGeneDir) && cd $(refGeneDir) && wget -nc $(genomeAssembly))


reference_genome: $(starBinDir)/STAR $(refGeneDir)/Homo_sapiens.GRCh38.94.gtf $(refGeneDir)/Homo_sapiens.GRCh38.dna.primary_assembly.fa
	mkdir -m 755 -p ${DNFA_starRoot}/STARIndex
	$(starBinDir)/STAR \
	  --runThreadN 12 \
	  --runMode genomeGenerate \
	  --genomeDir ${DNFA_starRoot}/STARIndex \
	  --genomeFastaFiles $(refGeneDir)/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	  --sjdbGTFfile $(refGeneDir)/Homo_sapiens.GRCh38.94.gtf \
	  --sjdbOverhang 75

