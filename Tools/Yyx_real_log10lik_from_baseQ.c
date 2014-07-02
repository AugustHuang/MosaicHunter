
/*********************************************************************
 * Copyright (C) 2013 Center for Bioinformatics, Peking University, China
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * General Public License for more details.
 *
 * This program is used to calculate L(r)=P(o|r,q)
 * using an iterative algorithm to traverse every possible 'real' alt count
 * Input includes ID(chr_pos), a column of ref baseQ, a column of alt baseQ
 *
 * Version: 0.5.1 (2014-01-14)\n");
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

#define PRED_TO_LOG10PROB(phred) (-1.0*(phred)/10.0)
#define POW10(x) (pow(10.0, x))

#define LOGZERO -1e100 

double log10_add_pow10(double a, double b);

double *yyx_calc_log10_prob_vec(const int*ref_quals, long ref_n, const int*alt_quals, long alt_n);



/**
 * @brief Computes log10(10^(a) + 10^(b))
 */
double log10_add_pow10(double a, double b)
{
/*  in R code:
	tmp_max <- max(vec)
	log10(sum(10^(vec-tmp_max))) + tmp_max
*/
//	printf("a = %f, b = %f, ", a, b);
	if (a > b) {
//		printf("ans = %f\n", a + log10(1+POW10(b-a)));
		return a + log10(1+POW10(b-a));
	} else {
//		printf("pow10 = %f\n", POW10(a-b));
//		printf("log10 = %f\n", log10(1+POW10(a-b)));
//		printf("ans = %f\n", b + log10(1+POW10(a-b)));
		return b + log10(1+POW10(a-b));
	}
}
/* end of log10_add_pow10() */

/**
 * @brief Computes log10 probability distribution (vector)
 *   of real alt count on the site
 *
 * using the algorithm as in R code:
	depth <- length(is_alt_vec)
	ans <- rep(-Inf, depth+1)
	ans[1] <- 0
	for(i in 1:depth){
		seq_wrong_log10_prob <- -baseQ_vec[i]/10
		seq_right_log10_prob <- log10(1 - 10^seq_wrong_log10_prob)
		alt_ref_prob <- c(seq_wrong_log10_prob, seq_right_log10_prob)
		if(is_alt_vec[i])
		{
			alt_ref_prob <- c(seq_right_log10_prob, seq_wrong_log10_prob)
		}
		for(j in (i+1):2)
		{
			ans[j] <- log10sum10exp(ans[c(j-1,j)] + alt_ref_prob)
		}
		ans[1] <- ans[1] + alt_ref_prob[2]
	}
	ans
 */
double *yyx_calc_log10_prob_vec(const int*ref_quals, long ref_n, const int*alt_quals, long alt_n)
{
	long depth = ref_n + alt_n;
	double *prob_vector = NULL;
	double log10_wrong, log10_right;
	long i, j;
	
	if((prob_vector = malloc((depth+2)*sizeof(double)))==NULL){
		fprintf(stderr, "Error: cannot allocate memory at %d\n", __LINE__);
		return NULL;
	}
	
	// initialization
	for(i=0; i<depth+2; i++){
		prob_vector[i] = LOGZERO;
	}
	// prob_vector[0] always be 0, 
	// and prob_vector[1:depth+1] is the probabilities for real_alt_count = 0:depth
	prob_vector[1] = 0.0;
	
	long now_depth = 0;
	// traverse every ref allele baseQ
	for(j=0; j<ref_n; j++){
		now_depth++;
//		printf("%d => %.5f\n", ref_quals[j], PRED_TO_LOG10PROB(ref_quals[j]));
		log10_wrong = PRED_TO_LOG10PROB(ref_quals[j]);
		log10_right = log10( 1 - POW10(log10_wrong) );
//		printf("log10 wrong = %.5f, right = %.5f\n", log10_wrong, log10_right);
		// as observing ref, wrong for real alt + 1, and right for real alt + 0
		for(i=now_depth+1; i>0; i--){
			prob_vector[i] = log10_add_pow10(prob_vector[i] + log10_right, prob_vector[i-1] + log10_wrong);
		}
	}
	
	// traverse every alt allele baseQ
	for(j=0; j<alt_n; j++){
		now_depth++;
		log10_wrong = PRED_TO_LOG10PROB(alt_quals[j]);
		log10_right = log10( 1 - POW10(log10_wrong) );
		// as observing alt, right for real alt + 1, and wrong for real alt + 0
		for(i=now_depth+1; i>0; i--){
			prob_vector[i] = log10_add_pow10(prob_vector[i] + log10_wrong, prob_vector[i-1] + log10_right);
		}
	}
	
	return prob_vector;
}
/* end of yyx_calc_log10_prob_vec() */


int main(int argc, char* argv[])
{
	char *ref_baseQ_str;
	char *alt_baseQ_str;
	int *ref_quals;
	int *alt_quals;
	long ref_n;
	long alt_n;
	
	char usage_str[1000];
	sprintf(usage_str, "\nUsage: cat <input> | %s <phred_shift> [max_depth] \n", argv[0]);
	strcat(usage_str, "\n");
	strcat(usage_str, " * This program is used to calculate L(r)=P(o|r,q)\n");
	strcat(usage_str, "   using an iterative algorithm to traverse every possible 'real' alt count\n");
	strcat(usage_str, " * Input includes ID(chr_pos), a column of ref baseQ, a column of alt baseQ\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  <phred_shift>: 33 for Phred+33, 64 for Phred+64\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  [max_depth]: \n");
	strcat(usage_str, "      max support depth (for memory allocate)\n");
	strcat(usage_str, "        [ default = 800000 ]\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  <input> format:\n");
	strcat(usage_str, "      field1: an ID string, such as chr:pos\n");
	strcat(usage_str, "      field2: a string of Phred baseQ for ref allele\n");
	strcat(usage_str, "      field3: a string of Phred baseQ for alt allele\n");
	strcat(usage_str, "        if no baseQ for ref or alt, you should fill the field with '|'\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "  output (stdout) format:\n");
	strcat(usage_str, "      ID +'\\t'+ depth +'\\t'+ log10 prob for 'real' alt count, separated by ','\n");
	strcat(usage_str, "\n");
	strcat(usage_str, "Version: 0.5.1 (2014-01-14)\n");
	strcat(usage_str, "Author: Adam Yongxin Ye @ CBI\n");
	strcat(usage_str, "\n");

	int phred_shift = 33;
	long max_support_depth = 800000;
	if(argc<2){
		fprintf(stderr, "%s", usage_str);
		return 1;
	}else{
//		phred_shift = argv[0]
		sscanf(argv[1], "%d", &phred_shift);
		if(argc>2){ sscanf(argv[2], "%ld", &max_support_depth); }
	}
	
	if((ref_baseQ_str = malloc((max_support_depth)*sizeof(char)))==NULL){
		fprintf(stderr, "Error: cannot allocate memory at %d\n", __LINE__);
		return 2;
	}
	if((alt_baseQ_str = malloc((max_support_depth)*sizeof(char)))==NULL){
		fprintf(stderr, "Error: cannot allocate memory at %d\n", __LINE__);
		free(ref_baseQ_str);
		return 2;
	}
	if((ref_quals = malloc((max_support_depth)*sizeof(int)))==NULL){
		fprintf(stderr, "Error: cannot allocate memory at %d\n", __LINE__);
		free(alt_baseQ_str);
		free(ref_baseQ_str);
		return 2;
	}
	if((alt_quals = malloc((max_support_depth)*sizeof(int)))==NULL){
		fprintf(stderr, "Error: cannot allocate memory at %d\n", __LINE__);
		free(ref_quals);
		free(alt_baseQ_str);
		free(ref_baseQ_str);
		return 2;
	}
	
	
	char chr_pos[1000];
	while(fscanf(stdin, "%s%s%s", chr_pos, ref_baseQ_str, alt_baseQ_str)!=EOF){
		long i = 0;
		// shift Phred for ref baseQ
		if(ref_baseQ_str[0]=='|'){
			ref_n = 0;
		}else{
			for(i=0; i<max_support_depth; i++){
				if(ref_baseQ_str[i]==0){
					break;
				}
				ref_quals[i] = (ref_baseQ_str[i]-phred_shift);
			}
			ref_n = i;
		}
		// shift Phred for alt baseQ
		if(alt_baseQ_str[0]=='|'){
			alt_n = 0;
		}else{
			for(i=0; i<max_support_depth; i++){
				if(alt_baseQ_str[i]==0){
					break;
				}
				alt_quals[i] = (alt_baseQ_str[i]-phred_shift);
			}
			alt_n = i;
		}
		
		// calculate the main purpose: log10 probability vector (length = depth + 2)
		double *log10_prob_vec = yyx_calc_log10_prob_vec(ref_quals, ref_n, alt_quals, alt_n);
		long depth = ref_n + alt_n;
		if(log10_prob_vec==NULL){
			fprintf(stdout, "0\t??? some error may have occured in yyx_calc_log10_prob_vec() ???\n");
		}else{
			fprintf(stdout, "%s\t%ld\t%.5f", chr_pos, depth, log10_prob_vec[1]);
			for(i=2; i<depth+2; i++){
				fprintf(stdout, ",%.5f", log10_prob_vec[i]);
			}
			fprintf(stdout, "\n");
			free(log10_prob_vec);
		}
	}
	
	free(alt_quals);
	free(ref_quals);
	free(alt_baseQ_str);
	free(ref_baseQ_str);
	
	return 0;
}
