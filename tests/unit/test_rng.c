#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "rng_dist.h"

#define N 100000
#define THRESHOLD 0.1

int test_gamma_sampling_ab(struct rng_state *rng)
{
	double a = 1000.0;
	double b = 5.0;
	
	double samples[N];
	double sum = 0;
	for(int i=0; i < N; i++) {
		double sample = gamma_sampling_ab(a, b, rng);
		samples[i] = sample;
		sum += sample;
	}
	
	double mean = sum / N;
	double variance = 0;
	for(int i=0; i < N; i++) {
		variance += pow(samples[i] - mean, 2);
	}
	variance /= N;

	assert(abs(mean - a/b) < THRESHOLD);
	assert(abs(variance - (a / pow(b, 2))) < THRESHOLD);

	return 0;
}

int test_beta_sampling(struct rng_state *rng)
{
	double a = 2.34;
	double b = 1.3;
	
	double samples[N];
	double sum = 0;
	for(int i=0; i < N; i++) {
		double sample = beta_sampling(a, b, rng);
		samples[i] = sample;
		sum += sample;
	}
	
	double mean = sum / N;
	double variance = 0;
	for(int i=0; i < N; i++) {
		variance += pow(samples[i] - mean, 2);
	}
	variance /= N;
	
	assert(abs(mean - a/(a+b)) < THRESHOLD);
	assert(abs(variance - (a*b / ((a+b)*(a+b)*(a+b+1)))) < THRESHOLD);

	return 0;
}

int test_binomial_sampling(struct rng_state *rng)
{
	int t = 100;
	double p = 0.1;
	
	int samples[N];
	int sum = 0;
	for(int i=0; i < N; i++) {
		int sample = binomial_sampling(t, p, rng);
		samples[i] = sample;
		sum += sample;
	}
	
	double mean = sum / N;
	double variance = 0;
	for(int i=0; i < N; i++) {
		variance += pow(samples[i] - mean, 2);
	}
	variance /= N;

	printf("%lf, %lf\n", mean, t*p);
	printf("%lf, %lf\n", variance, t*p*(1.0-p));
	
	assert(abs(mean - t*p) < 5);
	assert(abs(variance - t*p*(1.0-p)) < 10);

	return 0;
}

int main(void)
{
	unsigned int seed = time(NULL);
	struct rng_state *rng;
	rng = malloc(sizeof(struct rng_state));
	rng_init(rng, seed);

	test_gamma_sampling_ab(rng);
	printf("gamma_sampling_ab passed\n");
	test_beta_sampling(rng);
	printf("beta_sampling passed\n");
	test_binomial_sampling(rng);
	printf("binomial_sampling passed\n");
}


