/*********************************************************************
 *
 * * This program is just used to compute individual genotype posterior
 *   based on 2 files: individual genotype likelihood and prior,
 *   which is pretty simple, and to suit the recent framework
 * * <base_change_relative_rate> set the relative rate of base change
 *   relative to mosaic rate
 * * <genotype_log10lik> and <log10_prior> should be sorted
 *   on the first column ID(chr:pos:refBase:altBase)
 *
 * Version: 0.2.0 (2014-06-09)\n");
 * Author: Adam Yongxin Ye @ CBI\n");
 *
 *********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <unistd.h>

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

#define POW10(x) (pow(10.0, x))

#define LOGZERO -1e100 

#define LINE_LEN 1000
#define YYX_EOF_STR "EOF_Yyx_EOF_Yyx_EOF"

double mosaic_rate = 1e-7;
double minus_1_AF = 0.002;
double minus_2_AF = 1e-4;

char *chr_pos_delim;
char *chr_order[1000];
int chr_order_len = 0;
char individual_sex[5];

double base_change_matrix[16];   // row_idx * 4 + col_idx, row and col 0-A,1-C,2-G,3-T
double log10_add_pow10(double a, double b);

int compare_two_ID_chr_pos(const char *IDa, const char *IDb);

int base2idx(char base);

double *normalize_population_ACGT_allele_frequency(double *population_ACGT_allele_frequency, double minus_1_AF, double minus_2_AF);
double *individual_log10_prior(const char *chr_type, char individual_sex, const char *individual_ref_alt_base, const double *population_ACGT_allele_frequency, double mosaic_rate, double *tmp_4_double_vec);

double *yyx_calc_posterior(const double *log10_prior_vec, const double *log10_likelihood_vec, double *tmp_four_double_vec, double base_change_relative_rate);

char get_max_AF_base(const double *population_ACGT_allele_frequency, char base_excluded);
char get_higher_AF_base(const double *population_ACGT_allele_frequency, char base1, char base2);

int individual_genotyper_process_func(int num_str, const char **each_str);

int match_and_process_firstColumnSorted_files(int num_files, char **filenames, int (*compare_two_ID)(const char *IDa, const char *IDb), int (*process_func)(int num_str, const char **each_str), int process_error_should_stop );


/**
 * @brief Computes log10(10^(a) + 10^(b))
 */
double log10_add_pow10(double a, double b)
{
/*  in R code:
	tmp_max <- max(vec)
	log10(sum(10^(vec-tmp_max))) + tmp_max
*/
	if (a > b) {
		return a + log10(1+POW10(b-a));
//		return a + log1p(POW10(b-a))/log(10);
	} else {
		return b + log10(1+POW10(a-b));
//		return b + log1p(POW10(a-b))/log(10);
	}
}
/* end of log10_add_pow10() */



/**
 * @brief The compare function for sorting ID(chr_pos)
 *
 * compare ID, make YYX_EOF_STR the largest
 * if ID has '_', cut as chr_pos, str_cmp for chr and int_cmp for pos
 * otherwise, str_cmp ID
 */
int compare_two_ID_chr_pos(const char *IDa, const char *IDb)
{
	char modifiable_IDa[LINE_LEN];
	char modifiable_IDb[LINE_LEN];
	char *chr_a;
	char *chr_b;
	char *pos_str_a;
	char *pos_str_b;
	long pos_a, pos_b;
	int chr_idx_a, chr_idx_b;
	strcpy(modifiable_IDa, IDa);
	strcpy(modifiable_IDb, IDb);
	if(strcmp(IDa, YYX_EOF_STR)==0){
		if(strcmp(IDb, YYX_EOF_STR)==0){
			return 0;
		}else{   // IDb != YYX_EOF_STR
			return +1;   // make YYX_EOF_STR the largest
		}
	}else{   // IDa != YYX_EOF_STR
		if(strcmp(IDb, YYX_EOF_STR)==0){
			return -1;   // make YYX_EOF_STR the largest
		}else{   // IDb != YYX_EOF_STR
			if(strcmp(chr_pos_delim, "")!=0){
				pos_str_a = modifiable_IDa;
				chr_a = strsep(&pos_str_a, chr_pos_delim);
				sscanf(pos_str_a, "%ld", &pos_a);
				
				pos_str_b = modifiable_IDb;
				chr_b = strsep(&pos_str_b, chr_pos_delim);
				sscanf(pos_str_b, "%ld", &pos_b);
				
				if(strcmp(chr_a,chr_b)==0){
					return pos_a - pos_b;
				}else{
					for(chr_idx_a=0; chr_idx_a<chr_order_len; chr_idx_a++){
						if(strcmp(chr_a, chr_order[chr_idx_a])==0){
							break;
						}
					}
					for(chr_idx_b=0; chr_idx_b<chr_order_len; chr_idx_b++){
						if(strcmp(chr_b, chr_order[chr_idx_b])==0){
							break;
						}
					}
					return chr_idx_a - chr_idx_b;
				}
			}else{
				return strcmp(IDa, IDb);
			}
		}
	}
}
/* end of compare_two_ID_chr_pos() */



/**
 * @brief a test process_func()
 */
int test_process_func(int num_str, const char **each_str)
{
	int fi;
	for(fi=0; fi<num_str; fi++){
		printf("%d\t%s\n",  fi, each_str[fi]);
	}
	printf("\n");
	
	return 0;
}
/* end of test_process_func() */


/**
 * @brief array idx2base[] and function base2idx()
 */
const char *idx2base = "ACGT.";
int base2idx(char base){
	if(base=='A'){
		return 0;
	}else if(base=='C'){
		return 1;
	}else if(base=='G'){
		return 2;
	}else if(base=='T'){
		return 3;
	}else{
		return 4;
	}
}
/* end of base2idx() */


/**
 * @brief normalize 4 allele frequency, so that they sum to 1
 *
 * input population_ACGT_allele_frequency will be modified
 *
 * in fact, this is used to stablize mosaic_rate (1e-7)
 */
double *normalize_population_ACGT_allele_frequency(double *population_ACGT_allele_frequency, double minus_1_AF, double minus_2_AF)
{
//	double allele_freq_sum = 0;
	int bi;
	for(bi=0; bi<4; bi++){
		if(population_ACGT_allele_frequency[bi]>=0){
			// good, do nothing
		}else if(population_ACGT_allele_frequency[bi]==-1){
			population_ACGT_allele_frequency[bi] = minus_1_AF;
		}else if(population_ACGT_allele_frequency[bi]==-2){
			population_ACGT_allele_frequency[bi] = minus_2_AF;
		}else{
			fprintf(stderr, "Warning: unknown ref population allele frequency code '%f', I use '-2' instead\n", population_ACGT_allele_frequency[bi]);
			population_ACGT_allele_frequency[bi] = minus_2_AF;
		}
//		allele_freq_sum += population_ACGT_allele_frequency[bi];
	}
	
//	for(bi=0; bi<4; bi++){
//		 population_ACGT_allele_frequency[bi] /= allele_freq_sum;
//	}
	
	return population_ACGT_allele_frequency;
}
/* end of normalize_population_ACGT_allele_frequency() */

/**
 * @brief calculate individual's log10 prior based on population allele frequency
 *
 * you should first call normalize_population_ACGT_allele_frequency() to code '-1', '-2' and normalize
 *
 * output is log10 prior for ref-hom (or ref-hemi), het (or LOGZERO), alt-hom (or alt-hemi), mosaic
 */
double *individual_log10_prior(const char *chr_type, char individual_sex, const char *individual_ref_alt_base, const double *population_ACGT_allele_frequency, double mosaic_rate, double *tmp_4_double_vec)
{
	int individual_ref_base_idx = base2idx(individual_ref_alt_base[0]);
	int individual_alt_base_idx = base2idx(individual_ref_alt_base[1]);
	double ref_AF = population_ACGT_allele_frequency[individual_ref_base_idx];
	double alt_AF = population_ACGT_allele_frequency[individual_alt_base_idx];
	double log10_ref_AF;
	log10_ref_AF = log10(ref_AF/(ref_AF+alt_AF));   // normalize of the 2 allele/base
//	if(population_ACGT_allele_frequency[individual_ref_base_idx]>=0){
//		log10_ref_AF = log10(population_ACGT_allele_frequency[individual_ref_base_idx]);
//	}else if(population_ACGT_allele_frequency[individual_ref_base_idx]==-1){
//		log10_ref_AF = log10(minus_1_AF);
//	}else if(population_ACGT_allele_frequency[individual_ref_base_idx]==-2){
//		log10_ref_AF = log10(minus_2_AF);
//	}else{
//		fprintf(stderr, "Warning: unknown ref population allele frequency code '%f', I use '-2' instead\n", population_ACGT_allele_frequency[individual_ref_base_idx]);
//		log10_ref_AF = log10(minus_2_AF);
//	}
	double log10_alt_AF;
	log10_alt_AF = log10(alt_AF/(ref_AF+alt_AF));   // normalize of the 2 allele/base
//	if(population_ACGT_allele_frequency[individual_alt_base_idx]>=0){
//		log10_alt_AF = log10(population_ACGT_allele_frequency[individual_alt_base_idx]);
//	}else if(population_ACGT_allele_frequency[individual_alt_base_idx]==-1){
//		log10_alt_AF = log10(minus_1_AF);
//	}else if(population_ACGT_allele_frequency[individual_alt_base_idx]==-2){
//		log10_alt_AF = log10(minus_2_AF);
//	}else{
//		fprintf(stderr, "Warning: unknown alt population allele frequency code '%f', I use '-2' instead\n", population_ACGT_allele_frequency[individual_alt_base_idx]);
//		log10_alt_AF = log10(minus_2_AF);
//	}
	
	tmp_4_double_vec[3] = log10(mosaic_rate);
	double log10_not_mosiac = log10(1-mosaic_rate);
	if(chr_type[0]=='A' || (chr_type[0]=='X' && individual_sex=='F') ){   // autosome or (chrX and female
		tmp_4_double_vec[0] = 2*log10_ref_AF + log10_not_mosiac;
		tmp_4_double_vec[1] = log10(2)+log10_ref_AF+log10_alt_AF + log10_not_mosiac;
		tmp_4_double_vec[2] = 2*log10_alt_AF + log10_not_mosiac;
	}else if(chr_type[0]=='X' && individual_sex=='M'){   // chrX and male
		tmp_4_double_vec[0] = log10_ref_AF + log10_not_mosiac;
		tmp_4_double_vec[1] = LOGZERO;
		tmp_4_double_vec[2] = log10_alt_AF + log10_not_mosiac;
	}else if(chr_type[0]=='Y' && individual_sex=='M'){   // chrY and male
		tmp_4_double_vec[0] = log10_ref_AF + log10_not_mosiac;
		tmp_4_double_vec[1] = LOGZERO;
		tmp_4_double_vec[2] = log10_alt_AF + log10_not_mosiac;
	}else if(chr_type[0]=='Y' && individual_sex=='F'){   // chrY and female
		// for chrY and female, set all prior = 1
		// which genotyping results should be then filtered
		tmp_4_double_vec[0] = 0;
		tmp_4_double_vec[1] = 0;
		tmp_4_double_vec[2] = 0;
		tmp_4_double_vec[3] = 0;
	}else{
		fprintf(stderr, "Warning: unknown chr_type '%s' in individual_log10_prior()\n", chr_type);
	}
	
//	int i;
//	for(i=0; i<4; i++){
//		printf("[DEBUG] tmp_4_double_vec[%d] = %.5f\n", i, tmp_4_double_vec[i]);
//	}
	
	return tmp_4_double_vec;
}
/* end of individual_log10_prior() */

/**
 * @brief Computes genotype state posterior
 *
 * calculate posterior as R code below: 
 *
	ref_hom_posterior=ref_hom_prior*ref_hom_likelihood
	het_posterior=het_prior*het_likelihood
	alt_hom_posterior=alt_hom_prior*alt_hom_likelihood
	somatic_posterior=somatic_prior*somatic_likelihood
	sum_posterior=ref_hom_posterior+het_posterior+alt_hom_posterior+somatic_posterior
	c(ref_hom_posterior/sum_posterior,het_posterior/sum_posterior,alt_hom_posterior/sum_posterior,somatic_posterior/sum_posterior)
 */
double *yyx_calc_posterior(const double *log10_prior_vec, const double *log10_likelihood_vec, double *tmp_four_double_vec, double base_change_relative_rate)
{
	double alt_hom_log10_posterior = log10_prior_vec[2] + log10_likelihood_vec[2];
	double ref_hom_log10_posterior = log10_prior_vec[0] + log10_likelihood_vec[0];
	double log10_posterior_normalization_constant = log10_add_pow10(alt_hom_log10_posterior, ref_hom_log10_posterior);
	double het_log10_posterior = log10_prior_vec[1] + log10_likelihood_vec[1];
	log10_posterior_normalization_constant = log10_add_pow10(log10_posterior_normalization_constant, het_log10_posterior);
	double mosaic_log10_posterior = log10_prior_vec[3] + log10_likelihood_vec[3] + log10(base_change_relative_rate);
	log10_posterior_normalization_constant = log10_add_pow10(log10_posterior_normalization_constant, mosaic_log10_posterior);
	
	tmp_four_double_vec[0] = ref_hom_log10_posterior - log10_posterior_normalization_constant;
	tmp_four_double_vec[1] = het_log10_posterior - log10_posterior_normalization_constant;
	tmp_four_double_vec[2] = alt_hom_log10_posterior - log10_posterior_normalization_constant;
	tmp_four_double_vec[3] = mosaic_log10_posterior - log10_posterior_normalization_constant;
	
//	int i;
//	for(i=0; i<4; i++){
//		printf("log10_prior_vec[%d] = %.5f\n", i, log10_prior_vec[i]);
//	}
//	for(i=0; i<4; i++){
//		printf("log10_likelihood_vec[%d] = %.5f\n", i, log10_likelihood_vec[i]);
//	}
//	for(i=0; i<4; i++){
//		printf("tmp_four_double_vec[%d] = %.5f\n", i, tmp_four_double_vec[i]);
//	}
	
	return tmp_four_double_vec;
}
/* end of yyx_calc_posterior() */


/**
 * @brief find out base has highest pop_AF excluding the base_excluded
 */
char get_max_AF_base(const double *population_ACGT_allele_frequency, char base_excluded){
	double max_AF = -10;
	char max_AF_base = '.';
	int bi;
	for(bi=0; bi<4; bi++){
		if(idx2base[bi]!=base_excluded && population_ACGT_allele_frequency[bi]>max_AF){
			max_AF = population_ACGT_allele_frequency[bi];
			max_AF_base = idx2base[bi];
		}
	}
	
	return max_AF_base;
}
/* end of get_max_AF_base() */

/**
 * @brief find out which base has higher pop_AF
 *   (output base1 when equal)
 */
char get_higher_AF_base(const double *population_ACGT_allele_frequency, char base1, char base2){
	if(population_ACGT_allele_frequency[base2idx(base1)]>=population_ACGT_allele_frequency[base2idx(base2)]){
		return base1;
	}else{
		return base2;
	}
}
/* end of get_higher_AF_base() */


/**
 * @brief the process_func() for individual genotyper
 *
 * input must be 1(genotype_log10lik) + 1(prior_population_allele_frequency) = 2 strings
 *
 * read in data, construct CPD, calculate joint distribution (posterior) , and finally marginalize
 */
int individual_genotyper_process_func(int num_str, const char **each_str)
{
	if(num_str!=2){
		fprintf(stderr, "Error: the num_str for individual_genotyper_process_func() must be 2\n");
		return -1;
	}
//	test_process_func(num_str, each_str);
	
	char ID[1000] = "";
	char chr_pos[1000] = "";
	char *ID_ptr;
	char *chr;
	char *pos_str;
	char *refBase;
	char *altBase;
	double genotype_log10lik[4];
	char input_ref_alt_base[10];
	char ref_alt_base[10];
	char chr_type[10];
	double pop_ACGT_AF[5];
	// above four: [0] for ref-hom, [1] for het, [2] for alt-hom, [3] for mosaic
	int i;
	double posterior_log10[4];
	
	// read in each_str[0]: genotype_log10lik
	if(strcmp(each_str[0],"")==0){
		genotype_log10lik[0] = 1;
		genotype_log10lik[1] = 1;
		genotype_log10lik[2] = 1;
		genotype_log10lik[3] = 1;
		strcpy(input_ref_alt_base, "..");
	}else{
		sscanf(each_str[0], "%s\t%lf\t%lf\t%lf\t%lf", ID, &genotype_log10lik[0], &genotype_log10lik[1], &genotype_log10lik[2], &genotype_log10lik[3]);
		ID_ptr = ID;
//		printf("here1, ID = '%s'\n", ID);
		chr = strsep(&ID_ptr, chr_pos_delim);
		pos_str = strsep(&ID_ptr, chr_pos_delim);
//		printf("here2, chr = '%s', pos_str = '%s'\n", chr, pos_str);
		refBase = strsep(&ID_ptr, chr_pos_delim);
		altBase = strsep(&ID_ptr, chr_pos_delim);
//		printf("here3, refBase = '%s', altBase = '%s'\n", refBase, altBase);
		strcpy(chr_pos, chr);
		strcat(chr_pos, chr_pos_delim);
		strcat(chr_pos, pos_str);
		input_ref_alt_base[0] = refBase[0];
		input_ref_alt_base[1] = altBase[0];
		input_ref_alt_base[2] = '\0';
//		printf("here4, chr_pos = '%s', input_ref_alt_base = '%s'\n", chr_pos, input_ref_alt_base);
	}
	
	// read in each_str[1]: prior_population_allele_frequency
	if(strcmp(each_str[1],"")==0){
		if(strcmp(chr_pos, "")==0){
			fprintf(stderr, "Warning: no chr_pos info, and no valid ID for chr_pos, I just skipped it\n");
//			printf("here5\n");
			return 0;
		}
//		fprintf(stderr, "Warning: no prior info for ID = '%s'\n", chr_pos);
//		return 0;
//		strcpy(chr_pos, "");
		strcpy(chr_type, "A");
		pop_ACGT_AF[0] = -2;
		pop_ACGT_AF[1] = -2;
		pop_ACGT_AF[2] = -2;
		pop_ACGT_AF[3] = -2;
		pop_ACGT_AF[4] = -10;

		// assign prior's ref base to individual's ref
		if(input_ref_alt_base[0]!='.'){
			pop_ACGT_AF[base2idx(input_ref_alt_base[0])] = 1;
			fprintf(stderr, "Warning: no prior info for chr_pos = '%s', and I judge ref = individual's ref '%c'\n", chr_pos, input_ref_alt_base[0]);
		}else{
			fprintf(stderr, "Warning: no prior info and without any known base for chr_pos = '%s', I just skip it\n", chr_pos);
			return 0;
		}

	}else{
//		sscanf(each_str[1], "%s\t%lf\t%lf\t%lf\t%lf", chr_pos, &prior_log10[0], &prior_log10[1], &prior_log10[2], &prior_log10[3]);
		sscanf(each_str[1], "%s\t%s\t%lf\t%lf\t%lf\t%lf", chr_pos, chr_type, &pop_ACGT_AF[0], &pop_ACGT_AF[1], &pop_ACGT_AF[2], &pop_ACGT_AF[3]);
//		printf("here6, chr_type = '%s', chr_type[0] = '%c'\n", chr_type, chr_type[0]);
		if(chr_type[0]!='A' && chr_type[0]!='X' && chr_type[0]!='Y'){
			fprintf(stderr, "Warning: cannot recognize chr_type '%c', I just treat it 'A' (autosome)\n", chr_type[0]);
			strcpy(chr_type, "A");
		}
		pop_ACGT_AF[4] = -10;
	}
	
	// assign ref/alt base if necessary (troublesome)
	if(input_ref_alt_base[0]=='.'){   // no individual's ref/alt base
		ref_alt_base[0] = get_max_AF_base(pop_ACGT_AF, '.');
		ref_alt_base[1] = get_max_AF_base(pop_ACGT_AF, ref_alt_base[0]);
	}else if(input_ref_alt_base[1]=='.'){   // only individual's ref base, no alt base
		ref_alt_base[0] = input_ref_alt_base[0];
		ref_alt_base[1] = get_max_AF_base(pop_ACGT_AF, ref_alt_base[0]);
	}else{   // we know individual's ref/alt base
		ref_alt_base[0] = input_ref_alt_base[0];
		ref_alt_base[1] = input_ref_alt_base[1];
	}
	ref_alt_base[2] = '\0';
	
	
	// calculate log10 prior
	normalize_population_ACGT_allele_frequency(pop_ACGT_AF, minus_1_AF, minus_2_AF);
	double prior_log10[4];
	individual_log10_prior(chr_type, individual_sex[0], ref_alt_base, pop_ACGT_AF, mosaic_rate, prior_log10);
	
	// parse chr_pos, and select the corresponding base_change_relative_rate
	int row_idx = base2idx(ref_alt_base[0]);
	int col_idx = base2idx(ref_alt_base[1]);
	double base_change_relative_rate = 1;
	if(row_idx>3 || col_idx>3){
		fprintf(stderr, "Warning: refBase(%c) or altBase(%c) not in {A,C,G,T} for chr_pos = '%s', so I just set relative_rate=1\n", ref_alt_base[0], ref_alt_base[1], chr_pos);
	}else{
		base_change_relative_rate = base_change_matrix[row_idx*4+col_idx];
	}
	
	// compute posterior by multiplying prior and likelihood
	yyx_calc_posterior(prior_log10, genotype_log10lik, posterior_log10, base_change_relative_rate);
	
	
	// output format: ID(chr_pos) + chr_type + refBase/altBase + 4 genotype posterior
	printf("%s\t%s", chr_pos, chr_type);
	// ref/alt base
	printf("\t%c/%c", ref_alt_base[0], ref_alt_base[1]);
	if(strcmp(ref_alt_base, input_ref_alt_base)!=0){
		printf("(%c/%c)", input_ref_alt_base[0], input_ref_alt_base[1]);
	}
	if(individual_sex[0]=='F' && chr_type[0]=='Y'){
		printf("\t0\t0\t0\t0");
	}else{
		for(i=0; i<4; i++){   // genotype
			if(posterior_log10[i] < -1e8){
				printf("\t%.2e", posterior_log10[i]);
			}else{
				printf("\t%.5f", posterior_log10[i]);
			}
		}
	}
	printf("\n");
	
	return 0;
}
/* end of individual_genotyper_process_func() */



/**
 * @brief a wrapped framework to process matched line from sorted files
 *
 * open and read in each line from several sorted files, cut the first column as ID,
 * compare ID from each file (by compare_two_ID() ),
 * get the min_ID and put "" to those unmatched,
 * send these strings (may contain '\n') to process_func() for computing and outputing
 *
 * If process_func() returns not 0, show warning and proceed when process_error_should_stop=0
 * or stop, report error when process_error_should_stop=1
 */
int match_and_process_firstColumnSorted_files(int num_files, char **filenames, int (*compare_two_ID)(const char *IDa, const char *IDb), int (*process_func)(int num_str, const char **each_str), int process_error_should_stop )
{
	size_t malloc_line_len = LINE_LEN;
	const char *null_str = "";
	FILE **fileHandles;
	if((fileHandles = malloc((num_files)*sizeof(FILE *)))==NULL){
		fprintf(stderr, "Error: cannot allocate memory for fileHandles at %d\n", __LINE__);
		return 1;
	}
	int fi_0, fi_1;
	int has_error = 0;
	for(fi_0=0; fi_0<num_files; fi_0++){
		if((fileHandles[fi_0] = fopen(filenames[fi_0], "r"))==NULL){
			has_error = 1;
			fprintf(stderr, "Error: cannot open file '%s' for input\n", filenames[fi_0]);
//			break;
		}
	}
	int *should_read_next_line  = malloc((num_files)*sizeof(int));
	char **input_each_str = malloc((num_files)*sizeof(char*));
	int input_str_len;
	char **each_str = malloc((num_files)*sizeof(char*));
	char **each_ID_str = malloc((num_files)*sizeof(char*));
	if(should_read_next_line==NULL || each_str==NULL || input_each_str==NULL || each_ID_str==NULL){
		fprintf(stderr, "Error: cannot allocate memory for should_read_next_line, each_str, input_each_str or each_ID_str at %d\n", __LINE__);
		has_error = 1;
	}
	if(!has_error){
		for(fi_1=0; fi_1<num_files; fi_1++){
			should_read_next_line[fi_1] = 1;
			
			if((input_each_str[fi_1] = malloc((malloc_line_len)*sizeof(char)))==NULL){
				has_error = 1;
				fprintf(stderr, "Error: cannot allocate memory for input_each_str[%d]\n", fi_1);
//				break;
			}
			if((each_ID_str[fi_1] = malloc((malloc_line_len)*sizeof(char)))==NULL){
				has_error = 1;
				fprintf(stderr, "Error: cannot allocate memory for each_ID_str[%d]\n", fi_1);
//				break;
			}
		}
		
		char *min_ID_str;
		int not_all_eof = 1;
		int process_return_value;
		if(!has_error){
			while(not_all_eof){
				// read in new line
				for(fi_1=0; fi_1<num_files; fi_1++){
					if(should_read_next_line[fi_1]==1){
						if(strcmp(each_ID_str[fi_1], YYX_EOF_STR)==0){
							continue;
						}
						while((input_str_len = getline(&input_each_str[fi_1], &malloc_line_len, fileHandles[fi_1]))!=-1){
							if(input_str_len>0 && strcmp(input_each_str[fi_1], "")!=0 && strcmp(input_each_str[fi_1], "\n")!=0){
								break;
							}
						}
						if(input_str_len==-1){
							strcpy(each_ID_str[fi_1], YYX_EOF_STR);
						}else{
							sscanf(input_each_str[fi_1], "%s", each_ID_str[fi_1]);
						}
					}
				}
				
	//			printf("num_files = %d\n", num_files);
	//			for(fi_1=0; fi_1<num_files; fi_1++){
	//				printf("each_ID_str[%d] = %s\n", fi_1, each_ID_str[fi_1]);
	//			}
				// get min_ID_str
				min_ID_str = each_ID_str[0];
				for(fi_1=1; fi_1<num_files; fi_1++){
	//				printf("fi_1 = %d\n", fi_1);
					if(compare_two_ID(each_ID_str[fi_1],min_ID_str)<0){
						min_ID_str = each_ID_str[fi_1];
					}
				}
				if(strcmp(min_ID_str, YYX_EOF_STR)==0){
					not_all_eof = 0;
					continue;
				}
				
	//			printf("min_ID_str = %s\n", min_ID_str);
				// determine each_str and should_read_next_line
				for(fi_1=0; fi_1<num_files; fi_1++){
					if(compare_two_ID(each_ID_str[fi_1],min_ID_str)==0){
						each_str[fi_1] = input_each_str[fi_1];
						should_read_next_line[fi_1] = 1;
					}else{
						each_str[fi_1] = (char *)null_str;
						should_read_next_line[fi_1] = 0;
					}
				}
				
				// call process_func
				process_return_value = (*process_func)(num_files, (const char **)each_str);
				if(process_return_value!=0){
					if(process_error_should_stop){
						fprintf(stderr, "Error: process_func returns %d\n", process_return_value);
						break;
					}else{
						fprintf(stderr, "Warning: process_func returns %d\n", process_return_value);
					}
				}
			}

		}
		for(fi_0=num_files-1;fi_0>=0;fi_0--){
			if(each_ID_str[fi_0]!=NULL){
				free(each_ID_str[fi_0]);
			}
		}
	}
	free(should_read_next_line);
	free(each_str);
	free(input_each_str);
	free(each_ID_str);
	for(fi_0=num_files-1;fi_0>=0;fi_0--){
		if(fileHandles[fi_0]!=NULL){
			fclose(fileHandles[fi_0]);
		}
	}
	free(fileHandles);
	
	if(!has_error){
		return 0;
	}else{
		return -1;
	}
}
/* end of match_and_process_firstColumnSorted_files() */



int main(int argc, char* argv[])
{
	char usage_str[10000];
	sprintf(usage_str, "\nUsage: %s  <chr_pos_delim> <chr_order.list> \n", argv[0]);
	strcat(usage_str, "          <genotype_log10lik> <prior_allele_frequency>\n");
	strcat(usage_str, "          <individual_sex> [mosaic_rate] [minus_1_AF] [minus_2_AF] [base_change_relative_rate]\n");
	strcat(usage_str, "\n");
	strcat(usage_str, " * This program is just used to compute individual genotype posterior\n");
	strcat(usage_str, "   based on 2 files: individual genotype likelihood and prior,\n");
	strcat(usage_str, "   which is pretty simple, and to suit the recent framework\n");
	strcat(usage_str, " * <base_change_relative_rate> set the relative rate of base change\n");
	strcat(usage_str, "   relative to mosaic rate\n");
	strcat(usage_str, " * <genotype_log10lik> and <log10_prior> should be sorted\n");
	strcat(usage_str, "   on the first column ID(chr:pos:refBase:altBase)\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  <chr_pos_delim> : the char in ID separate chr and pos, eg. ':'\n");
	strcat(usage_str, "      use \"\" for just strcmp on ID\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  <chr_order.list> : a file contains chr order\n");
	strcat(usage_str, "      only consider the first column\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  <genotype_log10lik> format:\n");
	strcat(usage_str, "      ID(chr:pos:refBase:altBase) + log10 likelihood P(o|G)\n");
	strcat(usage_str, "        for genotype (ref-hom, het, alt-hom, mosaic)\n");
	strcat(usage_str, "        for each site on each line,\n");
	strcat(usage_str, "       can be generated by 'Yyx_genotype_log10lik_with_precalc_beta'\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  <prior_allele_frequency> format: (6 columns)\n");
	strcat(usage_str, "      ID(chr:pos) + A/X/Y (autosome/X/Y)\n");
	strcat(usage_str, "        + population allele frequency (A,C,G,T)\n");
	strcat(usage_str, "        for each site on each line,\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  <individual_sex>: individual's sex M/F (male/female)\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  [mosaic_rate]: mosaic rate\n");
	strcat(usage_str, "      the probability to see a site is a mosaic site in an individual\n");
	strcat(usage_str, "        [ default = 1e-7 ]\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  [minus_1_AF]: allele frequency for code '-1'\n");
	strcat(usage_str, "      the allele frequency for that in dbSNP but no MAF information\n");
	strcat(usage_str, "        [ default = 0.002 (1/500) ]\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  [minus_2_AF]: allele frequency for code '-2'\n");
	strcat(usage_str, "      the allele frequency for that not found in dbSNP\n");
	strcat(usage_str, "        [ default = 1e-4 ]\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  [base_change_relative_rate] format:\n");
	strcat(usage_str, "      4 * 4 table ('\\t' separated),\n");
	strcat(usage_str, "        the relative rate of change occurance,\n");
	strcat(usage_str, "        relative to the mosaic prior\n");
	strcat(usage_str, "        from row(A,C,G,T) to col(A,C,G,T)\n");
	strcat(usage_str, "      [ default: all = 1 ]\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  output (stdout) format: (7 columns, separated by '\\t')\n");
	strcat(usage_str, "      ID(chr:pos) + A/X/Y (autosome/X/Y) + refBase/altBase\n");
	strcat(usage_str, "      log10 posterior (ref-hom, het, alt-hom, mosaic)\n");
	strcat(usage_str, "        for each site on each line\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "Version: 0.2.0 (2014-06-09)\n");
	strcat(usage_str, "Author: Adam Yongxin Ye @ CBI\n");
	strcat(usage_str, "\n");
	
	char *chr_order_filename;
	char *base_change_filename = NULL;
	char **filenames;
	if(argc<5+1 || argc>9+1){
		fprintf(stderr, "%s", usage_str);
		return 1;
	}
	chr_pos_delim = argv[1];
	chr_pos_delim[1] = '\0';
	chr_order_filename = argv[2];
	filenames = argv+3;
	if(argc>5){
		sscanf(argv[5], "%s", individual_sex);
		if(strcmp(individual_sex,"M")!=0 && strcmp(individual_sex,"F")!=0){
			fprintf(stderr, "Error: cannot recognize individual's sex = '%c'\n", individual_sex[0]);
			return 2;
		}
	}
	if(argc>6){ sscanf(argv[6], "%lf", &mosaic_rate); }
	if(argc>7){ sscanf(argv[7], "%lf", &minus_1_AF); }
	if(argc>8){ sscanf(argv[8], "%lf", &minus_2_AF); }
	if(argc>9){ base_change_filename = argv[9];}
	
	int has_error = 0;
	size_t tmp_line_len = LINE_LEN;
	char *tmp_line = NULL;
	if((tmp_line = malloc((tmp_line_len)*sizeof(char)))==NULL){
		fprintf(stderr, "Error: cannot allocate memory for tmp_line\n");
		has_error = 1;
	}
	int return_value = -1;
	// read in base_change file
//	double base_change_matrix[16];   // row_idx * 4 + col_idx, row and col 0-A,1-C,2-G,3-T
	if(base_change_filename != NULL){
		FILE *base_change_fileHandle;
		char *ptr_for_strsep;
		if(!has_error){
			if((base_change_fileHandle = fopen(base_change_filename, "r"))==NULL){
				fprintf(stderr, "Error: cannot open file '%s' for input\n", base_change_filename);
				has_error = 1;
			}else{
				int row_idx = 0;
				int col_idx;
				while((getline(&tmp_line, &tmp_line_len, base_change_fileHandle))!=-1){
					if(strcmp(tmp_line, "")==0){
						continue;
					}
					ptr_for_strsep = tmp_line;
					for(col_idx=0; col_idx<4; col_idx++){
						base_change_matrix[row_idx*4+col_idx] = atof(ptr_for_strsep);
						strsep(&ptr_for_strsep, "\t");
					}
					row_idx++;
					if(row_idx>=4){
						break;
					}
				}
				fclose(base_change_fileHandle);
			}
		}
	}else{   // base_change_file not specified, so use default: all = 1
		int row_idx, col_idx;
		for(row_idx=0; row_idx<4; row_idx++){
			for(col_idx=0; col_idx<4; col_idx++){
				base_change_matrix[row_idx*4+col_idx] = 1;
			}
		}
	}
	// read in chr_order file
	if(!has_error){
		FILE *chr_order_fileHandle;
		if((chr_order_fileHandle = fopen(chr_order_filename, "r"))==NULL){
			fprintf(stderr, "Error: cannot open file '%s' for input\n", chr_order_filename);
			has_error = 1;
		}
		if(!has_error){
			while((getline(&tmp_line, &tmp_line_len, chr_order_fileHandle))!=-1){
				if(strcmp(tmp_line, "")==0){
					continue;
				}
				if((chr_order[chr_order_len] = malloc((tmp_line_len)*sizeof(char)))==NULL){
					fprintf(stderr, "Error: cannot allocate memory for chr_order[%d]\n", chr_order_len);
					has_error = 1;
					break;
				}
				sscanf(tmp_line, "%s", chr_order[chr_order_len]);
				chr_order_len++;
			}
			fclose(chr_order_fileHandle);
			
		//	printf("de_novo_rate = %.2e\n", de_novo_rate);
		//	printf("mosaic_rate = %.2e\n", mosaic_rate);
		//	return 0;
			
		//	return match_and_process_firstColumnSorted_files(6, filenames, compare_two_ID_chr_pos, test_process_func, 1 );
			if(!has_error){
				return_value = match_and_process_firstColumnSorted_files(2, filenames, compare_two_ID_chr_pos, individual_genotyper_process_func, 1 );
			}
			for(chr_order_len--; chr_order_len>=0; chr_order_len--){
				free(chr_order[chr_order_len]);
			}
		}else{
		}
	}else{
	}
	if(tmp_line != NULL){
		free(tmp_line);
		tmp_line = NULL;
	}
	return return_value;
}

