#!/usr/bin/make -s

TRINITY := /home/macmanes/trinityrnaseq_r2013-02-25
MUS := /media/macmanes/hd/flux/genomes/mus/Mus_musculus.GRCm38.71.cdna.all.fa
BCODES := /home/macmanes/Dropbox


RUN=run
READ1=left.fastq
READ2=right.fastq
TRIMMOMATIC ?= $(shell which 'trimmomatic-0.30.jar')

all: check trim correct assemble rsem

check:
	@echo "\n\n\n"###I am checking to see if you have all the dependancies installed.### "\n"
	command -v trimmomatic-0.30.jar >/dev/null 2>&1 || { echo >&2 "I require Trimmomatic but it's not installed.  Aborting."; exit 1; }
	@echo Trimmomatoc is Installed
	command -v SolexaQA.pl >/dev/null 2>&1 || { echo >&2 "I require SolexaQA.pl but it's not installed.  Aborting."; exit 1; }
	@echo SolexaQA is installed
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
		ILLUMINACLIP:$(BCODES)/barcodes.fa:2:40:15 \
		LEADING:$$TRIM TRAILING:$$TRIM SLIDINGWINDOW:4:$$TRIM MINLEN:25 ; \
		cat T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq > $(RUN)_left.$$TRIM.fastq ; \
		cat T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq > $(RUN)_right.$$TRIM.fastq ; \
		rm T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq ; \
	done

correct: hippo_left.$$TRIM.fastq hippo_right.$$TRIM.fastq
	@echo About to start error correction
	perl fastq-converter-v2.0.pl ./ ./ 1 #files MUST have fastq extension
	sed -i 's_[0-9]$_&/1_' hippo_left.fa #add /1 to ID reads as left
	sed -i 's_[0-9]$_&/2_' hippo_right.fa #add /2 to ID reads as right
	sed -i 's_^>.*[0-9]$_&/1_' hippo_left.q
	sed -i 's_^>.*[0-9]$_&/2_' hippo_right.q
	cat hippo_left.fa hippo_right.fa > both.fa
	cat hippo_left.q hippo_right.q > both.q
	$(REPTILE)/reptile-omp config.analy #Do error corection
	$(REPTILE)/utils/reptile_merger/reptile_merger both.fa both.reptile.err both.reptile.corr.fa #make error corrected fasta file
	grep -aA1 '/1' both.reptile.corr.fa > hippo.left.rept.corr.fa
	grep -aA1 '/2' both.reptile.corr.fa > hippo.right.rept.corr.fa

assemble:  hippo.left.rept.corr.fa hippo.right.rept.corr.fa
	$(TRINITY)/Trinity.pl --full_cleanup --SS_lib_type RF --min_kmer_cov 2 --seqType fa --JM 30G \
	--left hippo.left.rept.corr.fa --right hippo.right.rept.corr.fa --CPU 8 --output hippo
	
rsem: hippo.Trinity.fasta
	$(TRINITY)/util/RSEM_util/run_RSEM_align_n_estimate.pl --transcripts $< --seqType fq --left $(READ1) \
	--right $(READ2) --thread_count 8 --SS_lib_type RF -- --bowtie-chunkmbs 512


