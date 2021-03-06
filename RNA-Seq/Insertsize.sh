## following script is writen to calculate the insertion size of aligned pair-end reads in
## STAR generated BAM files.
## the following methods give very similar averge insert size statistics
## however bamtools does not provide STANDARD DEVIATION
## picard tool provides STANDARD DEVIATION but number is usually higher than 10%
# For mapping of DNA-seq data we can decide whether an alignment has insert size that's too big or too small by comparing it "expected" distribution of insert size.
# However, in RNA-seq the unsequenced portion of insert between the mates can contain a splice junctions.
# This means that we cannot simply calculate the insert size from an alignment to compare it with "expected" insert size.

########################################################
## Method 1: use bamtools to view stats from bam file ##
########################################################
## use Aligned.toTranscriptome.out.bam file as input
bamtools stats -in test6_S1_L003Aligned.toTranscriptome.out.bam -insert

# **********************************************
# Stats for BAM file(s):
# **********************************************

# Total reads:       33139892
# Mapped reads:      33139892	(100%)
# Forward strand:    16569946	(50%)
# Reverse strand:    16569946	(50%)
# Failed QC:         0	(0%)
# Duplicates:        0	(0%)
# Paired-end reads:  33139892	(100%)
# 'Proper-pairs':    33139892	(100%)
# Both pairs mapped: 33139892	(100%)
# Read 1:            16569946
# Read 2:            16569946
# Singletons:        0	(0%)
# Average insert size (absolute value): 206.381
# Median insert size (absolute value): 186


###################################################################
## Method 2: picard tool on bam files for insert size statistics ##
###################################################################
## install picard tool
sudo apt-get update
sudo apt-get install picard-tools
## find the pathway to picard tool #/home/suwu/.conda/pkgs/picard-2.18.16-0/share/picard-2.18.16-0/picard.jar
#/usr/share/java/picard.jar
## To test that you can run Picard tools, run the following command in your terminal application, providing either the full path to the picard.jar file
java -jar /usr/share/java/picard.jar -h
## analyze the bam file with picard tools for insert size statistics
## save the output files in a new folder picardinsert

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test1_S2_L001Aligned.toTranscriptome.out.bam \
      O=test1_1_insert_size_metrics.txt \
      H=test1_1_insert_size_histogram.pdf \
      M=0.5


java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test1_S2_L002Aligned.toTranscriptome.out.bam \
      O=test1_2_insert_size_metrics.txt \
      H=test1_2_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test1_S2_L003Aligned.toTranscriptome.out.bam \
      O=test1_3_insert_size_metrics.txt \
      H=test1_3_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test1_S2_L004Aligned.toTranscriptome.out.bam \
      O=test1_4_insert_size_metrics.txt \
      H=test1_4_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test2_S3_L001Aligned.toTranscriptome.out.bam \
      O=test2_1_insert_size_metrics.txt \
      H=test2_1_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test2_S3_L002Aligned.toTranscriptome.out.bam \
      O=test2_2_insert_size_metrics.txt \
      H=test2_2_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test2_S3_L003Aligned.toTranscriptome.out.bam \
      O=test2_3_insert_size_metrics.txt \
      H=test2_3_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test2_S3_L004Aligned.toTranscriptome.out.bam \
      O=test2_4_insert_size_metrics.txt \
      H=test2_4_insert_size_histogram.pdf \
      M=0.5

####
java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test3_S4_L001Aligned.toTranscriptome.out.bam \
      O=test3_1_insert_size_metrics.txt \
      H=test3_1_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test3_S4_L002Aligned.toTranscriptome.out.bam \
      O=test3_2_insert_size_metrics.txt \
      H=test3_2_insert_size_histogram.pdf \
      M=0.5
java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test3_S4_L003Aligned.toTranscriptome.out.bam \
      O=test3_3_insert_size_metrics.txt \
      H=test3_3_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test3_S4_L004Aligned.toTranscriptome.out.bam \
      O=test3_4_insert_size_metrics.txt \
      H=test3_4_insert_size_histogram.pdf \
      M=0.5

####
java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test4_S6_L001Aligned.toTranscriptome.out.bam \
      O=test4_1_insert_size_metrics.txt \
      H=test4_1_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test4_S6_L002Aligned.toTranscriptome.out.bam \
      O=test4_2_insert_size_metrics.txt \
      H=test4_2_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test4_S6_L003Aligned.toTranscriptome.out.bam \
      O=test4_3_insert_size_metrics.txt \
      H=test4_3_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test4_S6_L004Aligned.toTranscriptome.out.bam \
      O=test4_4_insert_size_metrics.txt \
      H=test4_4_insert_size_histogram.pdf \
      M=0.5

####
java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test5_S5_L001Aligned.toTranscriptome.out.bam \
      O=test5_1_insert_size_metrics.txt \
      H=test5_1_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test5_S5_L002Aligned.toTranscriptome.out.bam \
      O=test5_2_insert_size_metrics.txt \
      H=test5_2_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test5_S5_L003Aligned.toTranscriptome.out.bam \
      O=test5_3_insert_size_metrics.txt \
      H=test5_3_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test5_S5_L004Aligned.toTranscriptome.out.bam \
      O=test5_4_insert_size_metrics.txt \
      H=test5_4_insert_size_histogram.pdf \
      M=0.5

####
java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test6_S3_L001Aligned.toTranscriptome.out.bam \
      O=test6_1_insert_size_metrics.txt \
      H=test6_1_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test6_S3_L002Aligned.toTranscriptome.out.bam \
      O=test6_2_insert_size_metrics.txt \
      H=test6_2_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test6_S3_L003Aligned.toTranscriptome.out.bam \
      O=test6_3_insert_size_metrics.txt \
      H=test6_3_insert_size_histogram.pdf \
      M=0.5

java -jar /usr/share/java/picard.jar \
      CollectInsertSizeMetrics \
      I=/usr/local/DNFA-genfiles/data/r-extdata/test6_S3_L004Aligned.toTranscriptome.out.bam \
      O=test6_4_insert_size_metrics.txt \
      H=test6_4_insert_size_histogram.pdf \
      M=0.5

# Make a summary file from the tab-delimited summaries among all
# *metrics.txt files.  This will generate a tab-separated file
# that you can import into a spreadsheet program.  You will want
# to remove the text in cell 1A of the resulting spreadsheet.  It
# is spurious.
ls *metrics.txt | sort | xargs grep -A 1 MEDIAN_INSERT_SIZE | grep -v '\-\-' \
  | perl -e 'while(<>) { s/\.txt[:-]/\.txt\t/; print; } ' \
  > /tmp/table.tab && \
head -1 /tmp/table.tab > summary.tab && \
grep -v MEDIAN_INSERT_SIZE /tmp/table.tab | sort >> summary.tab

## check the statistics saved in output metics file
## head insert_size_metrics.txt
## htsjdk.samtools.metrics.StringHeader
# CollectInsertSizeMetrics HISTOGRAM_FILE=insert_size_histogram.pdf MINIMUM_PCT=0.5 INPUT=test6_S1_L003Aligned.toTranscriptome.out.bam OUTPUT=insert_size_metrics.txt    DEVIATIONS=10.0 METRIC_ACCUMULATION_LEVEL=[ALL_READS] INCLUDE_DUPLICATES=false ASSUME_SORTED=true STOP_AFTER=0 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
## htsjdk.samtools.metrics.StringHeader
# Started on: Thu Nov 15 20:59:10 EST 2018

## METRICS CLASS	picard.analysis.InsertSizeMetrics
## MEDIAN_INSERT_SIZE	MODE_INSERT_SIZE	MEDIAN_ABSOLUTE_DEVIATION	MIN_INSERT_SIZE	MAX_INSERT_SIZE	MEAN_INSERT_SIZE	STANDARD_DEVIATION	READ_PAIRS	PAIR_ORIENTATION	WIDTH_OF_10_PERCENT	WIDTH_OF_20_PERCENT	WIDTH_OF_30_PERCENT	WIDTH_OF_40_PERCENT	WIDTH_OF_50_PERCENT	WIDTH_OF_60_PERCENT	WIDTH_OF_70_PERCENT	WIDTH_OF_80_PERCENT	WIDTH_OF_90_PERCENT	WIDTH_OF_95_PERCENT	WIDTH_OF_99_PERCENT	SAMPLE	LIBRARY	READ_GROUP
## 189	154	44	35	15849	205.788198	78.000468	4132400	FR	17	35	53	71	89	109	135	169	251	349	595
