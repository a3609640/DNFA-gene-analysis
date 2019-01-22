##############################################
## Download SREBP1 ChIPSeq data from Encode ##
##############################################
cd ~/Documents/Bioinformatics_analysis/ChIP\ analysis/A549-SREBP1/Bam
# Download SREBF1 ChIPseq data in A549 cell line from Encode
# download SREBP1-A549-1.bam file
wget -c https://www.encodeproject.org/files/ENCFF224BFU/@@download/ENCFF224BFU.bam
# use md5sum for quality check after file download
md5sum A549-SREBP1-1.bam
# 963f1d263ab581902c394ccfebd1c293  A549-SREBP1-1.bam
ln -s ENCFF224BFU.bam A549-SREBP1-1.bam

wget -c https://www.encodeproject.org/files/ENCFF788IFN/@@download/ENCFF788IFN.bam
# quality check
md5sum ENCFF788IFN.bam
# 8c8124d093d66edcc9c319cd2defc9bb  ENCFF788IFN.bam
ln -s ENCFF788IFN.bam A549-SREBP1-2.bam

################################################################################
# Download ChIPSeq input data as control in A549 cell line from Encode
wget -c https://www.encodeproject.org/files/ENCFF189KJV/@@download/ENCFF189KJV.bam
# quality check
md5sum ENCFF189KJV.bam
# 64ac01c84b90eb1838ef558e5e504dbc  ENCFF189KJV.bam
ln -s ENCFF189KJV.bam A549-Input-1.bam

wget -c https://www.encodeproject.org/files/ENCFF419UIG/@@download/ENCFF419UIG.bam
md5sum ENCFF419UIG.bam
# 2c8a0d8e2ffc7ca65ade9f8fe9f0dfe1  ENCFF419UIG.bam
ln -s ENCFF419UIG.bam A549-Input-2.bam
# some error with this file when sorted: 'Input-A549-2.sorted.0015.bam' is truncated.

#  check the original file integrity
md5sum A549-SREBP1-1.bam A549-SREBP1-2.bam A549-Input-1.bam A549-Input-2.bam
# 963f1d263ab581902c394ccfebd1c293  A549-SREBP1-1.bam
# 8c8124d093d66edcc9c319cd2defc9bb  A549-SREBP1-2.bam
# 64ac01c84b90eb1838ef558e5e504dbc  A549-Input-1.bam
# 2c8a0d8e2ffc7ca65ade9f8fe9f0dfe1  A549-Input-2.bam

################################################################################
cd ~/Documents/Bioinformatics_analysis/ChIP\ analysis/MCF7-SREBP1/Bam
# Download SREBF1 chipseq data in MCF7 cell line from Encode (J. Michael Cherry, Stanford)
# download MCF7-SREBP1-1.bam file
wget -c https://www.encodeproject.org/files/ENCFF774XPS/@@download/ENCFF774XPS.bam
# quality check
md5sum ENCFF774XPS.bam
# e97bde1230c59824fad36d80e8ed9adc  ENCFF774XPS.bam
ln -s ENCFF774XPS.bam MCF7-SREBP1-1.bam

wget -c https://www.encodeproject.org/files/ENCFF008QCI/@@download/ENCFF008QCI.bam
# quality check
md5sum ENCFF008QCI.bam
# 6c2e59e7531eb2f6cf2855010f3cb22f ENCFF008QCI.bam
ln -s ENCFF008QCI.bam MCF7-SREBP1-2.bam

################################################################################
## Control Input pulldown file from the same cell line
wget -c https://www.encodeproject.org/files/ENCFF426RDP/@@download/ENCFF426RDP.bam
# quality check
md5sum ENCFF426RDP.bam
# b7b5004c4d1de084e602fdd2166f420d
ln -s ENCFF426RDP.bam MCF7-Input-1.bam

wget -c https://www.encodeproject.org/files/ENCFF291AFY/@@download/ENCFF291AFY.bam
# quality check
md5sum ENCFF291AFY.bam
# 84a7385de97ebc87cf22d1745fe03425
ln -s ENCFF291AFY.bam MCF7-Input-2.bam

md5sum MCF7-SREBP1-1.bam MCF7-SREBP1-2.bam MCF7-Input-1.bam MCF7-Input-2.bam

################################################################################
cd ~/Documents/Bioinformatics_analysis/ChIP\ analysis/K562-SREBP1/Bam
# Download SREBF1 chipseq data in K562 cell line from Encode
# download bam files
wget -c https://www.encodeproject.org/files/ENCFF405EME/@@download/ENCFF405EME.bam
# quality check
md5sum ENCFF405EME.bam
# d02ed64f85f40829fcedc97b5e649267  ENCFF405EME.bam
ln -s ENCFF405EME.bam K562-SREBP1-1.bam

wget -c https://www.encodeproject.org/files/ENCFF623PVZ/@@download/ENCFF623PVZ.bam
# quality check
md5sum ENCFF623PVZ.bam
# 90f97d5387748745914c1474a79e5a79  ENCFF623PVZ.bam
ln -s ENCFF623PVZ.bam K562-SREBP1-2.bam

################################################################################
# Download Control input pulldown file
wget -c https://www.encodeproject.org/files/ENCFF156FED/@@download/ENCFF156FED.bam
# quality check
md5sum ENCFF156FED.bam
# ee07408998c288c2b9b087589207caff  ENCFF156FED.bam
cp ENCFF156FED.bam K562-Input-1.bam

wget -c https://www.encodeproject.org/files/ENCFF577FNG/@@download/ENCFF577FNG.bam
# quality check
md5sum ENCFF577FNG.bam
# 8e1a12112216401ff4c7aaaa8314e6a1  ENCFF577FNG.bam
cp ENCFF577FNG.bam K562-Input-2.bam

md5sum K562-SREBP1-1.bam K562-SREBP1-2.bam K562-Input-1.bam K562-Input-2.bam
# d02ed64f85f40829fcedc97b5e649267  K562-SREBP1-1.bam
# 90f97d5387748745914c1474a79e5a79  K562-SREBP1-2.bam
# ee07408998c288c2b9b087589207caff  K562-Input-1.bam
# 8e1a12112216401ff4c7aaaa8314e6a1  K562-Input-2.bam


###############################################
## use samtools to sort downloaded bam files ##
###############################################
samtools sort A549-SREBP1-1.bam -o A549-SREBP1-1.sorted.bam
samtools sort A549-SREBP1-2.bam -o A549-SREBP1-2.sorted.bam
samtools sort A549-Input-1.bam -o A549-Input-1.sorted.bam
samtools sort A549-Input-2.bam -o A549-Input-2.sorted.bam

samtools sort MCF7-SREBP1-1.bam -o MCF7-SREBP1-1.sorted.bam
samtools sort MCF7-SREBP1-2.bam -o MCF7-SREBP1-2.sorted.bam
samtools sort MCF7-Input-1.bam -o MCF7-Input-1.sorted.bam
samtools sort MCF7-Input-2.bam -o MCF7-Input-2.sorted.bam

samtools sort K562-SREBP1-1.bam -o K562-SREBP1-1.sorted.bam
samtools sort K562-SREBP1-2.bam -o K562-SREBP1-2.sorted.bam
samtools sort K562-Input-1.bam -o K562-Input-1.sorted.bam
samtools sort K562-Input-2.bam -o K562-Input-2.sorted.bam

samtools index A549-SREBP1-1.sorted.bam
samtools index A549-SREBP1-2.sorted.bam
samtools index A549-Input-1.sorted.bam
samtools index A549-Input-2.sorted.bam

samtools index MCF7-SREBP1-1.sorted.bam
samtools index MCF7-SREBP1-2.sorted.bam
samtools index MCF7-Input-1.sorted.bam
samtools index MCF7-Input-2.sorted.bam

samtools index K562-SREBP1-1.sorted.bam
samtools index K562-SREBP1-2.sorted.bam
samtools index K562-Input-1.sorted.bam
samtools index K562-Input-2.sorted.bam

# merge bams file and add @RG and @CO
# Using samtools view -H <bamfile>, get @RG information for each of two bams.
samtools merge [-nr] [-h inh.sam] <out.bam> <in1.bam> <in2.bam> [...]
# merge replicate bam files for SREBP1 ChIP-Seq and Input respectively
samtools merge A549-merged-SREBP1.sorted.bam \
          A549-SREBP1-1.sorted.bam \
          A549-SREBP1-2.sorted.bam
samtools merge A549-merged-Input.sorted.bam \
          A549-Input-1.sorted.bam \
          A549-Input-2.sorted.bam

samtools merge MCF7-merged-SREBP1.sorted.bam \
          MCF7-SREBP1-1.sorted.bam \
          MCF7-SREBP1-2.sorted.bam
samtools merge MCF7-merged-Input.sorted.bam \
          MCF7-Input-1.sorted.bam \
          MCF7-Input-2.sorted.bam

samtools merge K562-merged-SREBP1.sorted.bam \
         K562-SREBP1-1.sorted.bam \
         K562-SREBP1-2.sorted.bam
samtools merge K562-merged-Input.sorted.bam \
         K562-Input-1.sorted.bam \
         K562-Input-2.sorted.bam

#########################################################
## MACS to call the peaks in SREBP1 ChIP-seq bam files ##
#########################################################
## install macs on linux machine
sudo apt-get install macs

# peak calling on merged bam files from merged SREBP1 sample and control sample
macs2 callpeak -t ./Bam/A549-merged-SREBP1.sorted.bam \
               -c ./Bam/A549-merged-Input.sorted.bam \
               -f BAM \
               -g hs \
               -n A549-merged-SREBP1-Input \
               --outdir MACS \
               -q 0.01 &> MACS/A549-merged-SREBP1-Input.log

macs2 callpeak -t ./Bam/MCF7-merged-SREBP1.sorted.bam \
               -c ./Bam/MCF7-merged-Input.sorted.bam \
               -f BAM \
               -g hs \
               -n MCF7-merged-SREBP1-Input \
               --outdir MACS \
               -q 0.01 &> MACS/MCF7-merged-SREBP1-Input.log

macs2 callpeak -t ./Bam/K562-SREBP1-1.sorted.bam \
               -c ./Bam/K562-Input-1.sorted.bam \
               -f BAM \
               -g hs \
               -n K562-SREBP1-1-Input-1 \
               --outdir MACS \
               -q 0.01 &> MACS/K562-SREBP1-1-Input-1.log
## the output files will be used for peak annotations with a R package

###################################################################
## generate bed files from ChIP-seq bam files for Motif analysis ##
###################################################################
## installation of bedtools
sudo apt-get install bedtools

cd ~/Documents/Bioinformatics_analysis/ChIP\ analysis/A549-SREBP1/Bed

bedtools bamtobed -i A549-SREBP1-1.bam > ../Bed/A549-SREBP1-1.bed
bedtools bamtobed -i A549-SREBP1-2.bam > ../Bed/A549-SREBP1-2.bed
bedtools bamtobed -i A549-Input-1.bam > ../Bed/A549-Input-1.bed
bedtools bamtobed -i A549-Input-2.bam > ../Bed/A549-Input-2.bed

bedtools bamtobed -i MCF7-SREBP1-1.bam > ../Bed/MCF7-SREBP1-1.bed
bedtools bamtobed -i MCF7-SREBP1-2.bam > ../Bed/MCF7-SREBP1-2.bed
bedtools bamtobed -i MCF7-Input-1.bam > ../Bed/MCF7-Input-1.bed
bedtools bamtobed -i MCF7-Input-2.bam > ../Bed/MCF7-Input-2.bed

bedtools bamtobed -i MCF7-SREBP1-1.bam > ../Bed/MCF7-SREBP1-1.bed
bedtools bamtobed -i MCF7-SREBP1-2.bam > ../Bed/MCF7-SREBP1-2.bed
bedtools bamtobed -i MCF7-Input-1.bam > ../Bed/MCF7-Input-1.bed
bedtools bamtobed -i MCF7-Input-2.bam > ../Bed/MCF7-Input-2.bed
# parallel < bamtobed.txt

##############################################################
## use homer to analyze the motifs in ChIP-seq files (.bed) ##
##############################################################
# install Homer
perl /home/suwu/Documents/Bioinformatics\ tools/Homer/configureHomer.pl -install homer
# Add the homer/bin directory to the executable path.  edit  ~/.bashrc file to include the line:
export PATH=$PATH:/home/suwu/Documents/Bioinformatics\ tools/Homer/bin/
# source /home/suwu/.bashrc

# To get a list of available packages:
perl /home/suwu/Documents/Bioinformatics\ tools/Homer/configureHomer.pl -list
perl /home/suwu/Documents/Bioinformatics\ tools/Homer/configureHomer.pl -update
# annotatePeaks.pl A549-merged-SREBP1-Input_peaks.narrowPeak hg38 > A549-SREBP1.txt

################################################################################
# Analyzing a ChIP-Seq experiment with one command
cd ~/Documents/Bioinformatics_analysis/ChIP\ analysis/A549-SREBP1/Bed
# Creating a "Tag Directory" with makeTagDirectory
makeTagDirectory ../Homer/TagDirectory/A549-Input A549-Input-1.bed A549-Input-2.bed  -format bed

analyzeChIP-Seq.pl ~/Documents/Bioinformatics_analysis/ChIP\ analysis/A549-SREBP1/Homer/TagDirectory/A549-SREBP1/ hg38 -i ~/Documents/Bioinformatics\ analysis/ChIP\ analysis/A549-SREBP1/Homer/TagDirectory/A549-Input -focus -A A549-SREBP1-1.bed A549-SREBP1-2.bed

sudo apt-get install PeakAnalyzer

################################################################################
cd ~/Documents/Bioinformatics_analysis/ChIP\ analysis/MCF7-SREBP1/MACS

makeTagDirectory ../Homer/TagDirectory/MCF7-Input MCF7-Input-1.bed MCF7-Input-2.bed  -format bed

analyzeChIP-Seq.pl ~/Documents/Bioinformatics_analysis/ChIP\ analysis/MCF7-SREBP1/Homer/TagDirectory/MCF7-SREBP1/ hg38 -i ~/Documents/Bioinformatics\ analysis/ChIP\ analysis/MCF7-SREBP1/Homer/TagDirectory/MCF7-Input -focus -A MCF7-SREBP1-1.bed MCF7-SREBP1-2.bed

################################################################################
cd ~/Documents/Bioinformatics\ analysis/ChIP\ analysis/K562-SREBP1/Bed

makeTagDirectory ../Homer/TagDirectory/K562-Input K562-Input-1.bed K562-Input-2.bed  -format bed

analyzeChIP-Seq.pl ~/Documents/Bioinformatics\ analysis/ChIP\ analysis/K562-SREBP1/Homer/TagDirectory/K562-SREBP1/ hg38 -i ~/Documents/Bioinformatics\ analysis/ChIP\ analysis/K562-SREBP1/Homer/TagDirectory/K562-Input -focus -A K562-SREBP1-1.bed K562-SREBP1-2.bed
