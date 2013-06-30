#!/usr/bin/make -s

SOLEXA := /home/macmanes/SolexaQA_v.2.1/SolexaQA_v.2.1
TRINITY := /home/macmanes/trinityrnaseq_r2013-02-25
TRIMMOMATIC := /home/macmanes/software
MUS := /media/macmanes/hd/flux/genomes/mus/Mus_musculus.GRCm38.71.cdna.all.fa
WD := /media/macmanes/raid/sequence-reads/broad_mrna/RT366294
BCODES := /home/macmanes/Dropbox
REPTILE := /home/macmanes/reptile-v1.1/reptile-v1.1

all: trim correct assemble

trim: $(WD)/raw_reads/tuco_hippocampus.1.fq $(WD)/raw_reads/tuco_hippocampus.2.fq
	@echo About to start trimming
	for TRIM in 10; do \
		java -Xmx30g -jar $(TRIMMOMATIC)/trimmomatic-0.30.jar PE -phred33 -threads 12 \
		$(WD)/raw_reads/tuco_hippocampus.1.fq \
		$(WD)/raw_reads/tuco_hippocampus.2.fq \
		T.$$TRIM.pp.1.fq \
		T.$$TRIM.up.1.fq \
		T.$$TRIM.pp.2.fq \
		T.$$TRIM.up.2.fq \
		ILLUMINACLIP:$(BCODES)/barcodes.fa:2:40:15 \
		LEADING:$$TRIM TRAILING:$$TRIM SLIDINGWINDOW:4:$$TRIM MINLEN:25 ; \
		cat T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq > hippo_left.$$TRIM.fastq ; \
		cat T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq > hippo_right.$$TRIM.fastq ; \
		rm T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq ; \
	done

correct:
	@echo About to start error correction
	#perl $(REPTILE)/utils/fastq-converter-v2.0.pl ./ ./ 1 #files MUST have fastq extension
	#sed -i 's_[0-9]$_&/1_' hippo_left.fa #add /1 to ID reads as left
	#sed -i 's_[0-9]$_&/2_' hippo_right.fa #add /2 to ID reads as right
	#sed -i 's_^>.*[0-9]$_&/1_' hippo_left.q
	#sed -i 's_^>.*[0-9]$_&/2_' hippo_right.q
	#cat hippo_left.fa hippo_right.fa > both.fa
	#cat hippo_left.q hippo_right.q > both.q
	$(REPTILE)/reptile-omp config.analy #Do error corection
	$(REPTILE)/reptile_merger both.fa both.reptile.err both.reptile.corr.fa #make error corrected fasta file
	grep -aA1 '/1' both.reptile.corr.fa > hippo.left.rept.corr.fa
	grep -aA1 '/2' both.reptile.corr.fa > hippo.right.rept.corr.fa

assemble:  
	$(TRINITY)/Trinity.pl --full_cleanup --SS_lib_type RF --min_kmer_cov 2 --seqType fa --JM 30G \
	--left hippo.left.rept.corr.fa --right hippo.right.rept.corr.fa --CPU 8 --output hippo
	
















