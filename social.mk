#!/usr/bin/make -s

SOLEXA := /home/macmanes/SolexaQA_v.2.1/SolexaQA_v.2.1
TRINITY := /home/macmanes/trinityrnaseq_r2013-02-25
TRIMMOMATIC := /home/macmanes/software
MUS := /media/macmanes/hd/flux/genomes/mus/Mus_musculus.GRCm38.71.cdna.all.fa
WD := /media/macmanes/raid/sequence-reads/broad_mrna/RT366294
BCODES := /home/macmanes/Dropbox
all: trim assemble

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
		cat T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq > hippo.left.$$TRIM.fq ; \
		cat T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq > hippo.right.$$TRIM.fq ; \
		rm T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq ; \
	done

assemble:  
	for TRIM in 10; do \
		$(TRINITY)/Trinity.pl --full_cleanup --SS_lib_type RF --min_kmer_cov 2 --seqType fq --JM 30G \
		--left hippo.left.$$TRIM.fq --right hippo.right.$$TRIM.fq --CPU 8 --output hippo ; \
	done
	

