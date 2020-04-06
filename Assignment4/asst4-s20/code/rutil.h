#ifndef RUTIL_H
/*
  C library implementing math and statistics functions for rat simulator.
*/

#include <stdlib.h>
#include <stdint.h>
#include <math.h>

/***************** Random number generation *******************/

/*
  Random number generator uses 64-bit arithmetic,
  but seed and random values guaranteed to fit in 32 bits
*/
typedef uint32_t random_t;

/* Standard parameters */
/* Default seed value */
#define DEFAULTSEED 618

/* Reinitialize seed based on list of seeds, where list has length len */
void reseed(random_t *seedp, random_t seed_list[], size_t len);

/* Generate double in range [0.0, upperlimit) */
double next_random_float(random_t *seedp, double upperlimit);

/* Select sample of size up to maxSample (without replacement) from list or size populationCount
 * Store result in array dest.  Array scratch must have space for sampleCount entries.
 * Returns minimum of populationCount and maxSample */
int sample(random_t *seedp, int *seq, int populationCount, int maxSample, int *dest, int *scratch);

/***************** Simulation                *******************/

/* Compute weight function */
double mweight(double val, double optval);

/* Compute imbalance between local and remote values */
/* Result < 0 when lcount > rcount and > 0 when lcount < rcount */
double imbalance(int lcount, int rcount);

/***************** Statistics                *******************/

/* Maximum of a set of elements */
double data_max(double *data, int n);

/* Sum of a set of elements */
double data_sum(double *data, int n);

/* Mean of a set of elements */
double data_mean(double *data, int n);

/* Standard deviation of a set of elements */
double data_stddev(double *data, int n);


/***************** Optimization                *******************/


/*
 * Given a set of values w_0, ..., w_n-1,
 * Find a partioning of these into k sets,
 * where each set has a contiguous range of values w_i, ..., w_i+j-1
 * Such if the cost of a partition is the sum of the values in the partition,
 * then the variance (and therefore the standard deviation) of these costs is minimized

 * The results are written to the array splits, where splits[m] is the
 * size of partition m
 */
void find_partition(int nweights, int npartitions, double *weights, int *splits);


#define RUTIL_H
#endif 
