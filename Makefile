
starSrcDir=$(DNFA_starRoot)/STAR-2.6.0a/source
starBinDir=$(DNFA_starRoot)/STAR-2.6.0a/bin/Linux_x86_64

${starBinDir}/STAR: ${starSrcDir}
	cd ${starSrcDir} && make 

${starSrcDir}: ${DNFA_starRoot}/2.6.0a.tar.gz
	$(shell cd $(DNFA_starRoot) && tar -xzf 2.6.0a.tar.gz)

$(DNFA_starRoot)/2.6.0a.tar.gz:
	$(shell cd $(DNFA_starRoot) && wget -nc https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz)

all: ${starBinDir}/STAR
