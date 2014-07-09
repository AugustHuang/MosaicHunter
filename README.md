<p>A script/tool for detecting postzygotic single-nucleotide mutations in human whole-genome sequencing data.</p>

<p>======Preparation</p>

<p>Make sure that you have installed the listed softwares, then added pre-installed softwares and the directory /your/MosaicHunter/directory/Tools in your PATH</p>

<p>In order to generate essential reference data and compile c, c++, java scripts, you should run this command once:</p>
<p>    seqpipe -m /your/MosaicHunter/directory/MosaicHunter.pipe preparation REFERENCE_DIR=/your/MosaicHunter/directory/Reference TOOLS_DIR=/your/MosaicHunter/directory/Tools</p>

<p>Pre-installed softwares required for the script:</p>
<p>    #SeqPipe: 0.4.12</p>
<p>    #BEDtools: 2.15.0</p>
<p>    #SAMtools: 0.1.18</p>
<p>    #FASTX-Toolkit: 0.0.13</p>
<p>    #Blat</p>
<p>    #fastasplitn</p>

<p>Reference data for the script: (Please put them into /your/MosaicHunter/directory/Reference)</p>
<p>    #human_g1k_v37.fasta (available at <a href="http://soms.nibs.ac.cn:6237/human_g1k_v37.fasta">http://soms.nibs.ac.cn:6237/human_g1k_v37.fasta</a>)</p>
<p>    #human_g1k_v37.genome</p>
<p>    #human_hg19.fasta (available at <a href="http://soms.nibs.ac.cn:6239/human_hg19.fasta">http://soms.nibs.ac.cn:6239/human_hg19.fasta</a>)</p>
<p>    #all_repeats.b37.bed</p>
<p>    #PAR.b37.bed</p>
<p>    #dbsnp_137.b37.SNP_AF.tsv (available at <a href="http://soms.nibs.ac.cn:6235/dbsnp_137.b37.SNP_AF.tsv">http://soms.nibs.ac.cn:6235/dbsnp_137.b37.SNP_AF.tsv</a>)</p>
<p>    #observed_in_common.bed</p>

<p>Tools for the script: (Please put them into /your/MosaicHunter/directory/Tools)</p>
<p>    #generate_beta_log10_val_file.r</p>
<p>    #count_homopolymer.cpp</p>
<p>    #myjoin</p>
<p>    #my_join.pl</p>
<p>    #PileupFilter.java</p>
<p>    #genotyper.pipe</p>
<p>    #Yyx_genotype_log10lik_with_precalc_beta.c</p>
<p>    #Yyx_real_log10lik_from_baseQ.c</p>
<p>    #Yyx_individual_genotyper.c</p>
<p>    #LoFreq_call.c</p>
<p>    #sam2fa.pl</p>
<p>    #blat_best.pipe</p>
<p>    #highest-score.pl</p>
<p>    #calculate-score-coverage-identity.pl</p>
<p>    #intersect_bed12.pipe</p>
<p>    #my.grep</p>
<p>    #trimBamByBlock.pl</p>
<p>    #strand_bias.R</p>
<p>    #allele_pos_dist.R</p>
<p>    #splitSamByAllele.pl</p>

<p>======Run</p>

<p>To identify pSNM sites from the whole-genome sequencing data, you can run this command: </p>
<p>    seqpipe -m /your/MosaicHunter/directory/MosaicHunter.pipe MosaicHunter REFERENCE_DIR=/your/MosaicHunter/directory/Reference TOOLS_DIR=/your/MosaicHunter/directory/Tools TEMP_DIR=/your/temp/directory INPUT_BAM=example.bam [INDEL_CNV_BED=example.bed] PROJECT_NAME=example GENDER=M THREAD_NUM=5</p>
<p>        [INPUT_BAM]: the path of your input .bam file, the .bam file should be sorted and indexed</p>
<p>        [INDEL_CNV_BED]: the path of a .bed file containing all the CNV and indel-flanking(+-5bp) regions which will be masked in our pipeline</p>
<p>        [PROJECT_NAME]: a string used as the prefix and suffix of the output files's name</p>
<p>        [GENDER]:  the gender of the subject, F or M</p>
<p>        [THREAD_NUM]: the maximum number of threads for running the script</p>

<p>Recommended pre-processing of the .bam file:</p>
<p>    1) Removing the duplicated, improper-paired, and multi-hit reads</p>
<p>    2) Removing the reads with more than three mismatches</p>
<p>    3) Processing the reads by GATK's indel realignment and base quality score recalibration</p>


<p>To change the running order and the parameters of the Bayesian genotyepr and the error filters, you can edit the the scripts of MosaicHunter in /MosaicHunter/MosaicHunter.pipe, according to the user manual of seqpipe.</p>

<p>======Output</p>

<p>The final list of the pSNM candidates could be found at MosaicHunter_[PROJECT_NAME]/[PROJECT_NAME].mosaic.final.tsv</p>
<p>    The colunms in the final list represent:</p>
<p>    1) chromosome</p>
<p>    2) position</p>
<p>    3) total depth</p>
<p>    4) reference nt</p>
<p>    5) alternative nt</p>
<p>    6) reference depth</p>
<p>    7) alternative depth</p>
<p>    8) -log10 of posterior probability of ref-hom genotype</p>
<p>    9) -log10 of posterior probability of het genotype</p>
<p>    10) -log10 of posterior probability of alt-hom genotype</p>
<p>    11) -log10 of posterior probability of mosaic genotype</p>
<p>    12) population allele fraction in dbSNP 137, -1 for annotated sites without information of allele fraction, -2 for unannotated sites</p>
<p>    13) sequence of +-500bp flanking regions</p>
