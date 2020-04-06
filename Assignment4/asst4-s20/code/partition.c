#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "crun.h"
#include <math.h>

#define THRESHOLD1 120.0
#define THRESHOLD2 300.0
#define DOUBLE_INF 1e30
/* This file is where you will write your partitioner */

/*
  This function should assign a zone id to every region in the graph.
  Zone IDs should range between 0 and Z-1.
  The updates should be made to the zone_id field of
  each entry in the region array

  Note that each region entry contains information about the position,
  node count, and edge count of the region.

  You should feel free to try out different partitioning schemes.
*/

int find_min_index(double *data, int size) {
  double min_data=DOUBLE_INF, min_index=-1;
  for (int i=0;i<size;i++) {
    if (data[i] < min_data) {
      min_data = data[i];
      min_index = i;
    }
  }
  return min_index;
}

void sort_by_key(double *weights, int *sorted_index, int size) {
  double *data = malloc(size*sizeof(double));
  memcpy(data, weights, size*sizeof(double));
  int min_index;
  for (int i=size-1;i>=0;i--) {
    min_index = find_min_index(data, size);
    data[min_index] = DOUBLE_INF;
    sorted_index[i] = min_index;
  }
}

void assign_zones(region_t *region_list, int nregion, int nzone) {
    // TODO.  This partitioner is very naive.  You can do better!

    // method3
    double *weights = malloc(nregion*sizeof(double));
    int *splits = malloc(nzone*sizeof(int));
    for (int rid=0;rid<nregion;rid++) {
      int nnode = region_list[rid].node_count;
      weights[rid] = nnode;
    }

    double std = data_stddev(weights,nregion);

    if (std <= THRESHOLD1) {
      find_partition(nregion,nzone,weights,splits);

      int rid=0;
      for (int zid=0;zid<nzone;zid++) {
        for (int s=0;s<splits[zid];s++) {
          region_list[rid++].zone_id=zid;
        }
      }
    } else if (std <= THRESHOLD2) {

      // greedy
      double *zone_weight = calloc(nzone, sizeof(double));
      int *sorted_index = malloc(nregion*sizeof(int));
      double weights_sum = data_sum(weights, nregion);
      for (int i=0;i<nregion;i++) {
        weights[i] = weights[i]/weights_sum;
      }
      sort_by_key(weights, sorted_index, nregion);

      int rid, min_index, last_index=0;
      for (int i=0;i<nregion;i++) {
        rid = sorted_index[i];
        min_index = find_min_index(zone_weight, nzone);
        // outmsg("rid: %d, minindex: %d, lastindex: %d, weights: %f, minweight: %f, lastweight: %f, diff: %f, ratio: %f\n", 
        //         rid, min_index, last_index, weights[rid], zone_weight[min_index], zone_weight[last_index], 
        //         zone_weight[last_index]-zone_weight[min_index], 
        //         (zone_weight[last_index]-zone_weight[min_index])/weights[rid]);

        if (zone_weight[last_index]-zone_weight[min_index]>weights[rid]) {
          zone_weight[min_index] += weights[rid];
          region_list[rid].zone_id = min_index;
          last_index = min_index;
        } else {
          zone_weight[last_index] += weights[rid];
          region_list[rid].zone_id = last_index;
        }
      }

    } else {
      // greedy
      double *zone_weight = calloc(nzone, sizeof(double));
      int *sorted_index = malloc(nregion*sizeof(int));
      double weights_sum = data_sum(weights, nregion);
      for (int i=0;i<nregion;i++) {
        weights[i] = weights[i]/weights_sum;
      }
      sort_by_key(weights, sorted_index, nregion);

      int rid, min_index;
      for (int i=0;i<nregion;i++) {
        rid = sorted_index[i];
        min_index = find_min_index(zone_weight, nzone);
        // outmsg("rid: %d, minindex: %d, lastindex: %d, weights: %f, minweight: %f, lastweight: %f, diff: %f, ratio: %f\n", 
        //         rid, min_index, last_index, weights[rid], zone_weight[min_index], zone_weight[last_index], 
        //         zone_weight[last_index]-zone_weight[min_index], 
        //         (zone_weight[last_index]-zone_weight[min_index])/weights[rid]);

        zone_weight[min_index] += weights[rid];
        region_list[rid].zone_id = min_index;

      }
    }



    // method 2
    // double *weights = malloc(nregion*sizeof(double));
    // int *splits = malloc(nzone*sizeof(int));
    // for (int rid=0;rid<nregion;rid++) {
    //   int nnode = region_list[rid].node_count;
    //   int nedge = region_list[rid].edge_count;
    //   weights[rid] = nnode;
    // }
    // find_partition(nregion,nzone,weights,splits);

    // int rid=0;
    // for (int zid=0;zid<nzone;zid++) {
    //   for (int s=0;s<splits[zid];s++) {
    //     region_list[rid++].zone_id=zid;
    //   }
    // }
    
    // method 1
    // for (int rid = 0; rid < nregion; rid++) {
	  //   region_list[rid].zone_id = rid % nzone;
    // }
}



