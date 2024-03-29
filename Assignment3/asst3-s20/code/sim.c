#include "crun.h"

#define THRESHOLD 8
#define TCOUNT omp_get_num_threads()
#define BLOCKW1 15
#define BLOCKH1 20
#define BLOCKW2 180
#define BLOCKH2 1
// #define SIZE1 8
#define SIZE2 3
#define SIZE3 64
// #define SIZE4 2400

#if DEBUG
/** USEFUL DEBUGGING CODE **/
static void show_weights(state_t *s) {
    int nid, eid;
    graph_t *g = s->g;
    int nnode = g->nnode;
    int *neighbor = g->neighbor;
    outmsg("Weights\n");
    for (nid = 0; nid < nnode; nid++) {
	int eid_start = g->neighbor_start[nid];
	int eid_end  = g->neighbor_start[nid+1];
	outmsg("%d: [sum = %.3f]", nid, compute_sum_weight(s, nid));
	for (eid = eid_start; eid < eid_end; eid++) {
	    outmsg(" %.3f", compute_weight(s, neighbor[eid]));
	}
	outmsg("\n");
    }
}
#endif

/*
  Compute all initial node counts according to rat population.
  Assumes that rat position array is initially zeroed out.
 */
static inline void take_census(state_t *s) {
    int *rat_position = s->rat_position;
    int *rat_count = s->rat_count;
    int nrat = s->nrat;
    int ri;
    for (ri = 0; ri < nrat; ri++) {
	rat_count[rat_position[ri]]++;
    }
}

/* Recompute all node weights */
static inline void compute_all_weights(state_t *s) {
    graph_t *g = s->g;
    double *node_weight = s->node_weight;
    START_ACTIVITY(ACTIVITY_WEIGHTS);
    
    // outmsg("tcount: %d\n", tcount);
    int tmpid;
    #pragma omp parallel
    {
        #pragma omp for nowait
        for (tmpid = 0; tmpid < g->bound; tmpid++){
            int nid = g->hubids[tmpid];
            int count = s->rat_count[nid];
            int outdegree = g->neighbor_start[nid+1] - g->neighbor_start[nid] - 1;
            int *start = &g->neighbor[g->neighbor_start[nid]+1];
            int i;
            double sum = 0.0;
            for (i = 0; i < outdegree; i++) {
                int lcount = s->rat_count[nid];
                int rcount = s->rat_count[start[i]];
                double r = imbalance(lcount, rcount);
                sum += r;
            }
            double ilf = BASE_ILF + ILF_VARIABILITY * (sum/outdegree);
            node_weight[nid] = mweight((double) count/s->load_factor, ilf);
        
        } 

        int width = g->width;
        int height = g->height;
        int blocksizew = BLOCKW1;
        int blocksizeh = BLOCKH1;
        int blockw = (width+blocksizew-1)/blocksizew;
        int blockh = (height+blocksizeh-1)/blocksizeh;
        int blocknum=blockw*blockh;

        #pragma omp for
        for (tmpid = 0; tmpid < blocknum; tmpid++){

            int blocki = tmpid % blockw;
            int blockj = tmpid / blockw;
            int tmpi, tmpj, tmph, tmpw;
            for (tmpj=0;tmpj<blocksizeh;tmpj++){
                for (tmpi=0;tmpi<blocksizew;tmpi++){
                    tmph = blockj*blocksizeh+tmpj;
                    tmpw = blocki*blocksizew+tmpi;
                    int nid = tmph*width+tmpw;
                    if (tmph < height && tmpw < width && !g->ishub[nid])
                    {
                    int count = s->rat_count[nid];
                    int outdegree = g->neighbor_start[nid+1] - g->neighbor_start[nid] - 1;
                    int *start = &g->neighbor[g->neighbor_start[nid]+1];
                    int i;
                    double sum = 0.0;
                    for (i = 0; i < outdegree; i++) {
                        int lcount = s->rat_count[nid];
                        int rcount = s->rat_count[start[i]];
                        double r = imbalance(lcount, rcount);
                        sum += r;
                    }
                    double ilf = BASE_ILF + ILF_VARIABILITY * (sum/outdegree);
                    node_weight[nid] = mweight((double) count/s->load_factor, ilf);
                    }
                }
            }
  
        } 
    }
   
    FINISH_ACTIVITY(ACTIVITY_WEIGHTS);
}

/* Precompute sums for each region */
static inline void find_all_sums(state_t *s) {
    graph_t *g = s->g;
    START_ACTIVITY(ACTIVITY_SUMS);

    int tmpid;

    #pragma omp parallel
    {   
        #pragma omp for nowait
        for (tmpid = 0; tmpid < g->bound; tmpid++){
            int nid = g->hubids[tmpid];
            double sum = 0.0;
            int eid;
            int stop = g->neighbor_start[nid+1];
            for (eid = g->neighbor_start[nid]; eid < stop; eid++) {
                sum += s->node_weight[g->neighbor[eid]];
                s->neighbor_accum_weight[eid] = sum;
            }
            s->sum_weight[nid] = sum;
        }

        int width = g->width;
        int height = g->height;
        int blocksizew = BLOCKW2;
        int blocksizeh = BLOCKH2;
        int blockw = (width+blocksizew-1)/blocksizew;
        int blockh = (height+blocksizeh-1)/blocksizeh;
        int blocknum=blockw*blockh;

        #pragma omp for schedule(dynamic, SIZE2)
        for (tmpid = 0; tmpid < blocknum; tmpid++){

            int blocki = tmpid % blockw;
            int blockj = tmpid / blockw;
            int tmpi, tmpj, tmph, tmpw;
            for (tmpj=0;tmpj<blocksizeh;tmpj++){
                for (tmpi=0;tmpi<blocksizew;tmpi++){
                    tmph = blockj*blocksizeh+tmpj;
                    tmpw = blocki*blocksizew+tmpi;
                    int nid = tmph*width+tmpw;
                    if (tmph < height && tmpw < width && !g->ishub[nid])
                    {
                    double sum = 0.0;
                    int eid;
                    int stop = g->neighbor_start[nid+1];
                    for (eid = g->neighbor_start[nid]; eid < stop; eid++) {
                        sum += s->node_weight[g->neighbor[eid]];
                        s->neighbor_accum_weight[eid] = sum;
                    }
                    s->sum_weight[nid] = sum;
                    }
                }
            }
        }

    }


    FINISH_ACTIVITY(ACTIVITY_SUMS);
}

/*
  Given list of increasing numbers, and target number,
  find index of first one where target is less than list value
*/

/*
  Linear search
 */
static inline int locate_value_linear(double target, double *list, int len) {
    int i;
    for (i = 0; i < len; i++)
	if (target < list[i])
	    return i;
    /* Shouldn't get here */
    return -1;
}
/*
  Binary search down to threshold, and then linear
 */
static inline int locate_value(double target, double *list, int len) {
    int left = 0;
    int right = len-1;
    while (left < right) {
	if (right-left+1 < BINARY_THRESHOLD)
	    return left + locate_value_linear(target, list+left, right-left+1);
	int mid = left + (right-left)/2;
	if (target < list[mid])
	    right = mid;
	else
	    left = mid+1;
    }
    return right;
}


/*
  This function assumes that node weights are already valid,
  and that have already computed sum of weights for each node,
  as well as cumulative weight for each neighbor
  Given list of integer counts, generate real-valued weights
  and use these to flip random coin returning value between 0 and len-1
*/
static inline int fast_next_random_move(state_t *s, int r) {
    int nid = s->rat_position[r];
    graph_t *g = s->g;
    random_t *seedp = &s->rat_seed[r];
    /* Guaranteed that have computed sum of weights */
    double tsum = s->sum_weight[nid];    
    double val = next_random_float(seedp, tsum);

    int estart = g->neighbor_start[nid];
    int elen = g->neighbor_start[nid+1] - estart;
    int offset = locate_value(val, &s->neighbor_accum_weight[estart], elen);
#if DEBUG
    if (offset < 0) {
	/* Shouldn't get here */
	outmsg("Internal error.  fast_next_random_move.  Didn't find valid move.  Target = %.2f/%.2f.\n",
	       val, tsum);
	return 0;
    }
#endif
    return g->neighbor[estart + offset];
}

static inline void compute_all_move(state_t *s, int bstart, int bcount){

    int ni, ri;
    graph_t *g = s->g;
    int nnode = g->nnode;
    START_ACTIVITY(ACTIVITY_NEXT);
    int* rat_position = s->rat_position;
    int* delta_rat_count = s->delta_rat_count;
    int* rat_count = s->rat_count;

    #pragma omp parallel for schedule(dynamic, SIZE3)
    for (ri = 0; ri < bcount; ri++) {
        int rid = ri+bstart;
        int onid = rat_position[rid];
        int nnid = fast_next_random_move(s, rid);
        rat_position[rid] = nnid;

        if (onid!=nnid){
            #pragma omp atomic
            delta_rat_count[onid] -= 1;
            #pragma omp atomic
            delta_rat_count[nnid] += 1;
        }
    }

    #pragma omp parallel for 
    for (ni = 0; ni < nnode; ni++) {
        rat_count[ni] += delta_rat_count[ni];
        // Clear count for future use
        delta_rat_count[ni] = 0;
    }

    FINISH_ACTIVITY(ACTIVITY_NEXT);
}


/* Process single batch */
static inline void do_batch(state_t *s, int bstart, int bcount) {
    find_all_sums(s);
    compute_all_move(s, bstart, bcount);
    /* Update weights */
    compute_all_weights(s);
}

static void batch_step(state_t *s) {
    int rid = 0;
    int bsize = s->batch_size;
    int nrat = s->nrat;
    int bcount;
    while (rid < nrat) {
	bcount = nrat - rid;
	if (bcount > bsize)
	    bcount = bsize;
	do_batch(s, rid, bcount);
	rid += bcount;
    }
}


double simulate(state_t *s, int count, update_t update_mode, int dinterval, bool display) {
    int i;
    /* Adjust bath size if not in bath mode */
    if (update_mode == UPDATE_SYNCHRONOUS)
	    s->batch_size = s->nrat;
    else if (update_mode == UPDATE_RAT)
	    s->batch_size = 1;

    /* Compute and show initial state */
    bool show_counts = true;
    double start = currentSeconds();
    take_census(s);
    compute_all_weights(s);
    if (display)
	    show(s, show_counts);

    for (i = 0; i < count; i++) {
        batch_step(s);
        if (display) {
            show_counts = (((i+1) % dinterval) == 0) || (i == count-1);
            show(s, show_counts);
        }
    }

    double delta = currentSeconds() - start;
    done();
    return delta;
}
