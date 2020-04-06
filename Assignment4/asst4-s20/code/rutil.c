/*
  C library implementing math and statistics functions for rat simulator.
*/

#include <stdio.h>
#include <stdbool.h>

#include "rutil.h"



/* Standard parameters */
#define GROUPSIZE 2147483647
#define MVAL  48271
#define VVAL  16807
#define INITSEED  418


static inline random_t rnext(random_t *seedp, random_t x) {
    uint64_t s = (uint64_t) *seedp;
    uint64_t xlong = (uint64_t) x;
    random_t val = ((xlong+1) * VVAL + s * MVAL) % GROUPSIZE;
    *seedp = (random_t) val;
    return val;
}

/* Reinitialize seed based on list of seeds, where list has length len */
void reseed(random_t *seedp, random_t seed_list[], size_t len) {
    *seedp = INITSEED;
    size_t i;
    for (i = 0; i < len; i++)
	rnext(seedp, seed_list[i]);
}

/* Generate double in range [0.0, upperlimit) */
double next_random_float(random_t *seedp, double upperlimit) {
    random_t val = rnext(seedp, 0);
    return ((double) val / (double) GROUPSIZE) * upperlimit;
}

/* Select sample of size up to maxSample (without replacement) from list or size populationCount */
/* Store result in array dest.  Array scratch must have space for sampleCount entries. */
/* Returns minimum of populationCount and maxSample */
int sample(random_t *seedp, int *seq, int populationCount, int maxSample, int *dest, int *scratch) {
    int i, idx, tmp;
    if (populationCount <= maxSample) {
	for (i = 0; i < populationCount; i++)
	    dest[i] = seq[i];
	return populationCount;
    }
    for (i = 0; i < maxSample; i++) {
	double w = next_random_float(seedp, 1.0);
	idx = i + (int) (w * (double) (populationCount-i));
	scratch[i] = idx;
	tmp = seq[idx];
	seq[idx] = seq[i];
	seq[i] = tmp;
	dest[i] = tmp;
    }
    for (i = maxSample-1; i >= 0; i--) {
	idx = scratch[i];
	tmp = seq[idx];
	seq[idx] = seq[i];
	seq[i] = tmp;
    }
    return maxSample;
}


/* Parameters for computing weights that guide next-move selection */
#define COEFF 0.4
#define OPTVAL 1.5

double mweight(double val, double optval) {
    double arg = 1.0 + COEFF * (val - optval);
    double lg = log(arg) * M_LOG2E;
    double denom = 1.0 + lg * lg;
    return 1.0/denom;
}

/* Compute imbalance between local and remote values */
/* Result < 0 when lcount > rcount and > 0 when lcount < rcount */
double imbalance(int lcount, int rcount) {
    if (lcount == 0 && rcount == 0)
	return 0.0;
    double sl = sqrt((double) lcount);
    double sr = sqrt((double) rcount);
    return (sr-sl)/(sr+sl);
}

/*** Statistics functions ***/


/* Maximum of a set of elements */
double data_max(double *data, int n) {
    double val = 0.0;
    int i;
    for (i = 0; i < n; ++i) {
	if (data[i] > val)
	    val = data[i];
    }
    return val;
}

/* Sum of a set of elements */
double data_sum(double *data, int n) {
    double sum = 0.0;
    int i;
    for (i = 0; i < n; ++i) {
        sum += data[i];
    }
    return sum;
}

/* Average of a set of elements */
double data_mean(double *data, int n) {
    if (n == 0)
	return 0.0;
    return data_sum(data, n) / (double) n;
}

/* Standard deviation of a set of elements */
/* Code adapted from https://www.programiz.com/c-programming/examples/standard-deviation */
double data_stddev(double *data, int n) {
    double mean2 = 0;
    int i;
    double mean = data_mean(data, n);
    for (i = 0; i < n; ++i) {
	double diff = data[i] - mean;
        mean2 += diff * diff;
    }
    return sqrt(mean2 / n);
}

/*** Optimization Functions ***/

/*
 * Below is the implementation of an efficient linear partitioner, based on dynamic programming.
 */


static double *lookupCost = NULL; // Cost of the optimal partition for (k, trimLength)
static int *lookupRlen = NULL;    // Size of the rightmost partition in optimal partition for (k, trimLength)

static double *allWeights = NULL;
static int allWeightCount = 0;

/* How many total entries were generated to table? */
static int entryCount = 0;

/* Helper functions */
static void setup_partition(int nweights, int npartitions, double *weights) {
    allWeights = weights;
    allWeightCount = nweights;
    lookupCost = calloc(nweights * npartitions+1, sizeof(double));
    lookupRlen = calloc(nweights * npartitions, sizeof(int));
    if (!lookupCost || !lookupRlen) {
	fprintf(stderr, "Couldn't allocate space for lookup tables\n");
	exit(1);
    }
    entryCount = 0;
}

static void finish_partition() {
    free(lookupCost);
    free(lookupRlen);
}

static inline int table_index(int k, int trimLength) {
    return (k-1) * allWeightCount + trimLength;
}
    
static inline bool checkTable(int k, int trimLength) {
    return lookupRlen[table_index(k, trimLength)] != 0;
}

/* What is the cost of a subrange of weights */
static inline double segmentCost(int leftIndex, int length) {
    double sum = 0.0;
    int i;
    for (i = leftIndex; i < leftIndex + length; i++)
	sum += allWeights[i];
    return sum * sum;
}

/* At end, construct optimal partition */
static void construct_splits(int npartitions, int *splits) {
    int trimLength = 0;
    int k = npartitions;
    while (k > 0) {
	if (!checkTable(k, trimLength)) {
	    fprintf(stderr, "Partitioner failed.  Couldn't find table entry (%d, %d)\n", k, trimLength);
	    exit(0);
	}
	int rlen = lookupRlen[table_index(k, trimLength)];
	splits[k-1] = rlen;
	trimLength += rlen;
	k --;
    }
}

/*
 * Construct table of optimal subsolutions
 * Each indexed by (k', trimLength), where k' is the size of a partition
 * and trimLength indicates how much to reduce the upper part of the weights
 */
static void buildTable(int k, int trimLength) {
    if (checkTable(k, trimLength))
	return;
    
    int n = allWeightCount - trimLength;
    int idx = table_index(k, trimLength);
    int bestRlen;
    double bestCost;
    
    if (k == 1) {
	/* First n nodes should be put into single partition */
	bestRlen = n;
	bestCost = segmentCost(0, n);
    } else {
	bestRlen = 0;
	bestCost = -1.0;
	int rlen;
	// Look at ways to create new partition on right 
	for (rlen = 1; rlen <= n-k+1; rlen++) {
	    double segCost = segmentCost(n-rlen, rlen);
	    buildTable(k-1, trimLength+rlen);
	    int lidx = table_index(k-1, trimLength+rlen);
	    double restCost = lookupCost[lidx];
	    double cost = restCost + segCost;
	    if (bestCost < 0 || bestCost > cost) {
		bestCost = cost;
		bestRlen = rlen;
	    }
	}
    }
    lookupCost[idx] = bestCost;
    lookupRlen[idx] = bestRlen;
    entryCount++;
}

/* Here's the code for the partitioner */
void find_partition(int nweights, int npartitions, double *weights, int *splits) {
    /* Deal with trivial cases */
    if (npartitions == 1) {
	splits[0] = nweights;
	return;
    }
    if (npartitions >= nweights) {
	int i;
	for (i = 0; i < npartitions; i++) {
	    splits[i] = i < nweights ? 1 : 0;
	}
	return;
    }
    setup_partition(nweights, npartitions, weights);
    buildTable(npartitions, 0);
    construct_splits(npartitions, splits);
    finish_partition();
}
