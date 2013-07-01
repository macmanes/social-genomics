#!/usr/bin/make -rRsf

###########################################
###        -usage 'social.mk READ1=/location/or/read1.fastq READ2=/location/of/read2.fastq'
###
###         -Make sure your Trinity base directory 
###         	is set properly
###         -Make sure barcode file is located with
###           BCODES= tag
############################################
TRINITY := /home/macmanes/trinityrnaseq_r2013-02-25
BCODES := /home/macmanes/Dropbox/barcodes.fa
CONFIG:= /home/macmanes/Dropbox/config.analy



##### No Editing should be necessary below this line  #####






RUN=run
READ1=left.fastq
READ2=right.fastq
TRIMMOMATIC ?= $(shell which 'trimmomatic-0.30.jar')

all: check trim correct merge assemble rsem

check:
	@echo "\n\n\n"###I am checking to see if you have all the dependancies installed.### "\n"
	command -v trimmomatic-0.30.jar >/dev/null 2>&1 || { echo >&2 "I require Trimmomatic but it's not installed.  Aborting."; exit 1; }
	@echo Trimmomatic is Installed
	command -v $(TRINITY/Trinity.pl) >/dev/null 2>&1 || { echo >&2 "I require Trinity but it's not installed.  Aborting."; exit 1; }
	@echo Trinity is Installed
	command -v fastq-converter-v2.0.pl >/dev/null 2>&1 || { echo >&2 "I require fastq-converter-v2.0.pl (Reptile package) but it's not installed.  Aborting."; exit 1; }
	command -v reptile-omp >/dev/null 2>&1 || { echo >&2 "I require reptile-omp but it's not installed.  Aborting."; exit 1; }
	command -v reptile_merger >/dev/null 2>&1 || { echo >&2 "I require reptile_merger but it's not installed.  Aborting."; exit 1; }
	@echo Reptile is installed"\n"

trim: $(READ1) $(READ2)
	@echo About to start trimming
	for TRIM in 10; do \
		java -Xmx30g -jar $(TRIMMOMATIC) PE -phred33 -threads 12 \
		$(READ1) \
		$(READ2) \
		T.$$TRIM.pp.1.fq \
		T.$$TRIM.up.1.fq \
		T.$$TRIM.pp.2.fq \
		T.$$TRIM.up.2.fq \
		ILLUMINACLIP:$(BCODES):2:40:15 \
		LEADING:$$TRIM TRAILING:$$TRIM SLIDINGWINDOW:4:$$TRIM MINLEN:25 ; \
		cat T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq > $(RUN)_left.$$TRIM.fastq ; \
		cat T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq > $(RUN)_right.$$TRIM.fastq ; \
		rm T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq ; \
	done

correct: $(RUN)_left.$$TRIM.fastq $(RUN)_right.$$TRIM.fastq
	@echo About to start error correction
	perl fastq-converter-v2.0.pl ./ ./ 1 #files MUST have fastq extension
	sed -i 's_[0-9]$_&/1_' $(RUN)_left.fa #add /1 to ID reads as left
	sed -i 's_[0-9]$_&/2_' $(RUN)_right.fa #add /2 to ID reads as right
	sed -i 's_^>.*[0-9]$_&/1_' $(RUN)_left.q
	sed -i 's_^>.*[0-9]$_&/2_' $(RUN)_right.q
	cat $(RUN)_left.fa $(RUN)_right.fa > both.fa
	cat $(RUN)_left.q $(RUN)_right.q > both.q
	reptile-omp $(CONFIG) #Do error corection

merge: both.reptile.err
	reptile_merger both.fa $< both.reptile.corr.fa #make error corrected fasta file
	grep -aA1 '/1' both.reptile.corr.fa > $(RUN).left.rept.corr.fa
	grep -aA1 '/2' both.reptile.corr.fa > $(RUN).right.rept.corr.fa
	

assemble:  $(RUN).left.rept.corr.fa $(RUN).right.rept.corr.fa
	$(TRINITY)/Trinity.pl --full_cleanup --SS_lib_type RF --min_kmer_cov 2 --seqType fa --JM 30G \
	--left $(RUN).left.rept.corr.fa --right $(RUN).right.rept.corr.fa --CPU 8 --output $(RUN)
	
rsem: $(RUN).Trinity.fasta
	$(TRINITY)/util/RSEM_util/run_RSEM_align_n_estimate.pl --transcripts $< --seqType fq --left $(READ1) \
	--right $(READ2) --thread_count 8 --SS_lib_type RF -- --bowtie-chunkmbs 512


