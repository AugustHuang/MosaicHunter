MosaicHunter
============

A script/tool for detecting postzygotic single-nucleotide mutations in human whole-genome sequencing data.

======Preparation

Make sure that you have installed the listed softwares, then added pre-installed softwares and the directory /your/MosaicHunter/directory/Tools in your PATH

In order to generate essential reference data and compile c, c++, java scripts, you should run this command once:
	seqpipe -m /your/MosaicHunter/directory/MosaicHunter.pipe preparation REFERENCE_DIR=/your/MosaicHunter/directory/Reference TOOLS_DIR=/your/MosaicHunter/directory/Tools

Pre-installed softwares required for the script:
	#bedtools: 2.15.0
	#samtools: 0.1.18
	#fastx_toolkit: 0.0.13
	#seqpipe: 0.4.12
	#blat
	#fastasplitn

Reference data for the script: (Please put them into /your/MosaicHunter/directory/Reference)
	#human_g1k_v37.fasta (available at http://soms.nibs.ac.cn:6237/human_g1k_v37.fasta)
	#human_g1k_v37.genome
	#human_hg19.fasta (available at http://soms.nibs.ac.cn:6239/human_hg19.fasta)
	#all_repeats.b37.bed
	#PAR.b37.bed
	#dbsnp_137.b37.SNP_AF.tsv (available at http://soms.nibs.ac.cn:6235/dbsnp_137.b37.SNP_AF.tsv)
	#observed_in_common.bed

Tools for the script: (Please put them into /your/MosaicHunter/directory/Tools)
	#generate_beta_log10_val_file.r
	#count_homopolymer.cpp
	#myjoin
	#my_join.pl
	#PileupFilter.java
	#genotyper.pipe
	#Yyx_genotype_log10lik_with_precalc_beta.c
	#Yyx_real_log10lik_from_baseQ.c
	#Yyx_individual_genotyper.c
	#LoFreq_call.c
	#sam2fa.pl
	#blat_best.pipe
	#highest-score.pl
	#calculate-score-coverage-identity.pl
	#intersect_bed12.pipe
	#my.grep
	#trimBamByBlock.pl
	#strand_bias.R
	#allele_pos_dist.R
	#splitSamByAllele.pl
	
======Run

To identify pSNM sites from the whole-genome sequencing data, you can run this command: 
	seqpipe -m /your/MosaicHunter/directory/MosaicHunter.pipe MosaicHunter REFERENCE_DIR=/your/MosaicHunter/directory/Reference TOOLS_DIR=/your/MosaicHunter/directory/Tools TEMP_DIR=/your/temp/directory INPUT_BAM=example.bam INDEL_CNV_BED=example.bed PROJECT_NAME=example GENDER=M THREAD_NUM=5
		[INPUT_BAM]: the path of your input .bam file, the .bam file should be sorted and indexed.		 
		[INDEL_CNV_BED]: the path of a .bed file containing all the CNV and indel-flanking(+-5bp) regions which will be masked in our pipeline.  
		[PROJECT_NAME]: a string used as the prefix and suffix of the output files's name
		[GENDER]:  the gender of the subject, F or M
		[THREAD_NUM]: the maximum number of threads for running the script
		
	Recommended pre-processing of the .bam file: 
	1) Removing the duplicated, improper-paired, and multi-hit reads
	2) Removing the reads with more than three mismatches
	2) Processing the reads by GATK's indel realignment and base quality score recalibration
	
To change the running order and the parameters of the Bayesian genotyepr and the error filters, you can edit the the scripts of MosaicHunter in /MosaicHunter/MosaicHunter.pipe, according to the user manual of seqpipe.

======Output

The final list of the pSNM candidates could be found at MosaicHunter_[PROJECT_NAME]/[PROJECT_NAME].mosaic.final.tsv
	The colunms in the final list represent: 1) chromosome
											 2) position
											 3) total depth
											 4) reference nt
											 5) alternative nt
											 6) reference depth
											 7) alternative depth
											 8) -log10 of posterior probability of ref-hom genotype
											 9) -log10 of posterior probability of het genotype
											 10) -log10 of posterior probability of alt-hom genotype
											 11) -log10 of posterior probability of mosaic genotype
											 12) population allele fraction in dbSNP 137, -1 for annotated sites without information of allele fraction, -2 for unannotated sites
											 13) sequence of +-500bp flanking regions
