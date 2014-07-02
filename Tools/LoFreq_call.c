/*********************************************************************
 *
 * Copyright (C) 2011, 2012 Genome Institute of Singapore
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
 *********************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <unistd.h>

#define PHREDQUAL_TO_PROB(phred) (pow(10.0, -1.0*(phred)/10.0))
#define PHREDQUAL_VALID_RANGE(phred) ((phred)>1 && (phred)<100)
#define PROB_TO_PHREDQUAL(prob) ((-10.0 * log10(prob)))

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

#define LOGZERO -1e100 
#define FLOAT_DELTA 1e-32


double log_sum(double log_a, double log_b);
double log_diff(double log_a, double log_b);
double probvec_tailsum(double *probvec, int tail_startindex,
                       int probvec_len);
double *naive_calc_prob_dist(const int *quals, int N, int K);
double *pruned_calc_prob_dist(const int *quals, int N, int K);


/**
 * @brief Computes log(exp(log_a) + exp(log_b))
 *
 * Taken from util.h of FAST source code:
 * http://www.cs.cornell.edu/~keich/FAST/fast.tar.gz
 * and using log1p
 */
double
log_sum(double log_a, double log_b)
{
    if (log_a > log_b) {
        return log_a + log1p(exp(log_b-log_a));
    } else {
        return log_b + log1p(exp(log_a-log_b));
    }
}
/* end of log_sum() */


/**
 * @brief Computes log(exp(log_a) - exp(log_b))
 *
 * Adapted from log_sum above and scala/breeze/numerics logDiff
 * See also http://stackoverflow.com/questions/778047/we-know-log-add-but-how-to-do-log-subtract
 *
 */
double
log_diff(double log_a, double log_b)
{
    if (log_a >= log_b) {
        return log_a + log1p(- exp(log_b-log_a));
    } else {
        return log_b + log1p(- exp(log_a-log_b));
    }
}
/* end of log_diff() */


/**
 * @brief Computes sum of probvec values (log space) starting from (including)
 * tail_startindex to (excluding) probvec_len
 *
 */
double
probvec_tailsum(double *probvec, int tail_startindex, int probvec_len)
{
    double tailsum;
    int i;

    tailsum = probvec[tail_startindex];
    for (i=tail_startindex+1; i<probvec_len; i++) {
        tailsum = log_sum(tailsum, probvec[i]);
    }

    return tailsum;
}
/* end of probvec_tailsum() */


/**
 * @brief Computes the likelihood of tail_index (log space) in probvec
 *
 */
double
probvec_likelihood(double *probvec, int tail_index)
{
    return probvec[tail_index];
}
/* end of probvec_likelihood() */


/**
 * @brief FIXME:missing-doc
 *
 *
 */
double *
naive_calc_prob_dist(const int *quals, int N, int K)
{
    double *probvec = NULL;
    double *probvec_prev = NULL;
    double *probvec_swp = NULL;

    int n;

    if (NULL == (probvec = malloc((N+1) * sizeof(double)))) {
        fprintf(stderr, "FATAL: couldn't allocate memory at %d\n",
                __LINE__);
        return NULL;
    }
    if (NULL == (probvec_prev = malloc((N+1) * sizeof(double)))) {
        fprintf(stderr, "FATAL: couldn't allocate memory at %d\n",
                __LINE__);
        return NULL;
    }

    /* init */
    probvec_prev[0] = 0.0; /* 0.0 = log(1.0) */

    for (n=1; n<N+1; n++) {
        int k;
        double pn = PHREDQUAL_TO_PROB(quals[n-1]);
        double log_pn = log(pn);
        double log_1_pn = log1p(-pn);

        assert(PHREDQUAL_VALID_RANGE(quals[n-1]));

        k = 0;
        probvec[k] = probvec_prev[k] + log_1_pn;

        for (k=1; k<n; k++) {
            probvec[k] = log_sum(probvec_prev[k] + log_1_pn,
                                 probvec_prev[k-1] + log_pn);
        }
        k = n;
        probvec[k] = probvec_prev[k-1] + log_pn;
		
		/* swap */
		if (n < N) {
			probvec_swp = probvec;
			probvec = probvec_prev;
			probvec_prev = probvec_swp;
		}
    }
	
    free(probvec_prev);    
    return probvec;
}
/* end of naive_prob_dist */


/**
 * @brief FIXME:missing-doc
 *
 *
 */
double *
pruned_calc_prob_dist(const int *quals, int N, int K)
{
    double *probvec = NULL;
    double *probvec_prev = NULL;
    double *probvec_swp = NULL;
    int n;

    if (NULL == (probvec = malloc((K+1) * sizeof(double)))) {
        fprintf(stderr, "FATAL: couldn't allocate memory at %d\n",
                __LINE__);
        return NULL;
    }
    if (NULL == (probvec_prev = malloc((K+1) * sizeof(double)))) {
        fprintf(stderr, "FATAL: couldn't allocate memory at %d\n",
                __LINE__);
        return NULL;
    }

    /* init */
    probvec_prev[0] = 0.0; /* log(1.0) */

    for (n=1; n<=N; n++) {
        int k;
        double pvalue;
        double pn = PHREDQUAL_TO_PROB(quals[n-1]);
        double log_pn = log(pn);
        double log_1_pn = log1p(-pn); /* 0.0 = log(1.0) */
        
        /* test for valid phred quality boundaries */
        assert(PHREDQUAL_VALID_RANGE(quals[n-1]));

        if(n < K) {
            probvec_prev[n] = LOGZERO;
        }

        for (k=MIN(n,K-1); k>=1; k--) {
            assert(probvec_prev[k]<=0.0 && probvec_prev[k-1]<=0.0);
            probvec[k] = log_sum(probvec_prev[k] + log_1_pn,
                                 probvec_prev[k-1] + log_pn);            
        }
        k = 0;
        assert(probvec_prev[k]<=0.0);
        probvec[k] = probvec_prev[k] + log_1_pn;

        if (n==K) {
            probvec[K] = probvec_prev[K-1] + log_pn;
            /* FIXME check here as well */

        } else if (n > K) { 
            assert(probvec_prev[K]<=0.0 && probvec_prev[K-1]<=0.0);
            probvec[K] = log_sum(probvec_prev[K], probvec_prev[K-1]+log_pn);
            pvalue = exp(probvec[K]);     
		}
		
		/* swap */
		if (n < N) {
			probvec_swp = probvec;
			probvec = probvec_prev;
			probvec_prev = probvec_swp;
		}
    }

 free_and_exit:
    free(probvec_prev);    
    return probvec;
}
/* end of pruned_calc_prob_dist */


/**
 * @brief
 * 
 * Call per pileup column
 *
 * If pvalue was not computed (always insignificant) its value
 * will be set to DBL_MAX
 * 
 */
void
snpcaller(const int *phred_quals, const int num_phred_quals, 
          const int num_alt_phred_quals, const int naive,
		  double *snp_pvalue, double *snp_likelihood)
{
    double *probvec = NULL;
	*snp_pvalue = DBL_MAX;
	*snp_likelihood = DBL_MAX;

	if (naive) {
	    probvec = naive_calc_prob_dist(phred_quals, num_phred_quals,
									   num_alt_phred_quals);
		*snp_pvalue = exp(probvec_tailsum(probvec, num_alt_phred_quals, num_phred_quals+1));
/* 
		*snp_likelihood = exp(probvec_likelihood(probvec, num_alt_phred_quals));
*/
	}
	else {
		probvec = pruned_calc_prob_dist(phred_quals, num_phred_quals,
										num_alt_phred_quals);
		*snp_pvalue = exp(probvec_tailsum(probvec, num_alt_phred_quals, num_alt_phred_quals+1));
/*
		*snp_likelihood = exp(probvec_likelihood(probvec, num_alt_phred_quals));
*/
	}

/*	
	int i=0;
	for (i=0; i<=num_phred_quals; i++) {
		if (naive) {
			fprintf(stderr, "DEBUG(%s:%s():%d): pvalue=%g for num_alt_phred_quals %d\n", 
					__FILE__, __FUNCTION__, __LINE__, 
					exp(probvec_tailsum(probvec, i, num_phred_quals+1)), i);
		}
		else {
			fprintf(stderr, "DEBUG(%s:%s():%d): pvalue=%g for num_alt_phred_quals %d\n", 
					__FILE__, __FUNCTION__, __LINE__, 
					exp(probvec_tailsum(probvec, i, i+1)), i);
		}
	}
*/
	
	if (NULL != probvec) {
        free(probvec);
    }
}
/* end of snpcaller() */


/**
 * @brief Computes likelihoods of three genotypes (ref_hom, alt_hom, het)
 * based on sequencing qualities
 *
 */
void
genotype_likelihood(const int *phred_quals, const int num_phred_quals,
		   const int num_alt_phred_quals, double *ref_hom,
		   double *alt_hom, double *het)
{
	int i;
	double pn;
	double log_pn;
	double log_1_pn;
	
	*ref_hom=0;
	*alt_hom=0;
	*het=0;
	
	for (i=0; i<num_phred_quals-num_alt_phred_quals; i++)
	{
		pn = PHREDQUAL_TO_PROB(phred_quals[i]);
        log_pn = log(pn);
        log_1_pn = log1p(-pn);
		*ref_hom = *ref_hom + log_1_pn;
		*alt_hom = *alt_hom + log_pn;
		*het = *het + log(0.5);
	}
	
	for (i=num_phred_quals-num_alt_phred_quals; i<num_phred_quals; i++)
	{
		pn = PHREDQUAL_TO_PROB(phred_quals[i]);
        log_pn = log(pn);
        log_1_pn = log1p(-pn);
		*ref_hom = *ref_hom + log_pn;
		*alt_hom = *alt_hom + log_1_pn;
		*het = *het + log(0.5);
	}
}
/* end of genotype_likelihood() */


int main(int argc, char* argv[])
{
	char line[102400];
	int quals[102400];

	int opt;
	int i;
	int naive;
	int qbase;
	
	int num_phred_quals;
    int num_alt_phred_quals;
	
	double pvalue;
	double likelihood;
	double ref_hom;
	double alt_hom;
	double het;
	
	if (argc!=5)
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: cat <input> | %s -n <0/1> -q <int>\n", argv[0]);
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "    -n <0/1>      (Required)    naive (but slower)\n");
		fprintf(stderr, "    -q <int>      (Required)    ascii_base\n");
		fprintf(stderr, "<input> format:\n");
		fprintf(stderr, "    field1: a string of ascii-based phred scores\n");
		fprintf(stderr, "    field2: the number of reference and alternative bases\n");
		fprintf(stderr, "    field3: the number of alternative bases\n");
		fprintf(stderr, "<output> format:\n");
		fprintf(stderr, "    field1: the phred score of probability in LoFreq (Wilm et al.)\n");
		fprintf(stderr, "    field2: the natural logarithm of likelihood for genotype ref_hom\n");
		fprintf(stderr, "    field3: the natural logarithm of likelihood for genotype het\n");
		fprintf(stderr, "    field4: the natural logarithm of likelihood for genotype alt_hom\n\n");
		return 1;
	}
	
	while((opt=getopt(argc, argv, "n:q:"))!=-1)
	{
		if (opt=='n')
		{
			sscanf(optarg, "%d", &naive);
		}
		if (opt=='q')
		{
			sscanf(optarg, "%d", &qbase);
		}
	}
	
	while(fscanf(stdin, "%s%d%d", line, &num_phred_quals, &num_alt_phred_quals)!=EOF)
	{
		for (i=0; i<num_phred_quals; i++)
		{
			quals[i]=line[i]-qbase;
		}
		snpcaller(quals, num_phred_quals, num_alt_phred_quals, naive, &pvalue, &likelihood);
		genotype_likelihood(quals, num_phred_quals, num_alt_phred_quals, &ref_hom, &alt_hom, &het);
		if (naive) {
			fprintf(stdout, "%f\t%f\t%f\t%f\n", fabs(PROB_TO_PHREDQUAL(pvalue)), ref_hom, het, alt_hom);
		}
		else {
			fprintf(stdout, "%f\t%f\t%f\t%f\n", fabs(PROB_TO_PHREDQUAL(pvalue)), ref_hom, het, alt_hom);
		}
	}	
	
	return 0;
}