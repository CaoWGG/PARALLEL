
#include "crun.h"


/* Compute ideal load factor (ILF) for node */
static inline double neighbor_ilf(state_t *s, int nid) {
    graph_t *g = s->g;
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
    double ilf = BASE_ILF + 0.5 * (sum/outdegree);
    return ilf;
}

/* Compute weight for node nid */
static inline double compute_weight(state_t *s, int nid) {
    int count = s->rat_count[nid];
    double ilf = neighbor_ilf(s, nid);
    return mweight((double) count/s->load_factor, ilf);
}


/* Recompute all node counts according to rat population */
/*
  Function only called at start of simulation, at which point
  every zone has complete rat information.  Can therefore
  have every zone update every node count.
*/
static inline void take_census(state_t *s) {
    graph_t *g = s->g;
    int nnode = g->nnode;
    int *rat_position = s->rat_position;
    int *rat_count = s->rat_count;
    int nrat = s->nrat;

    memset(rat_count, 0, nnode * sizeof(int));
    int ri;
    for (ri = 0; ri < nrat; ri++) {
	rat_count[rat_position[ri]]++;
    }

}

/* Recompute all node weights */
static inline void compute_all_weights(state_t *s) {
    int nid;
    graph_t *g = s->g;
    double *node_weight = s->node_weight;
    START_ACTIVITY(ACTIVITY_WEIGHTS);
    for (nid = 0; nid < g->nnode; nid++)
	node_weight[nid] = compute_weight(s, nid);
    FINISH_ACTIVITY(ACTIVITY_WEIGHTS);

}

static inline void compute_cur_weights(state_t *s) {
    int nid;
    graph_t *g = s->g;
    double *node_weight = s->node_weight;
    START_ACTIVITY(ACTIVITY_WEIGHTS);

    for (int tmpid=0;tmpid<g->local_node_count;tmpid++) {
        nid = g->local_node_list[tmpid];
	    node_weight[nid] = compute_weight(s, nid);
    }

    FINISH_ACTIVITY(ACTIVITY_WEIGHTS);

}

/* In synchronous or batch mode, can precompute sums for each region */
static inline void find_all_sums(state_t *s) {
    graph_t *g = s->g;
    START_ACTIVITY(ACTIVITY_SUMS);
    // TODO: It doesn't make sense to compute the weights for nodes that are not
    // in the local zone
    int nid, eid;
#if MPI
    for (int tmpid=0;tmpid<g->local_node_count;tmpid++) {
    nid = g->local_node_list[tmpid];
#else
    for (nid=0;nid<g->nnode;nid++){
#endif
    double sum = 0.0;
	for (eid = g->neighbor_start[nid]; eid < g->neighbor_start[nid+1]; eid++) {
	    sum += s->node_weight[g->neighbor[eid]];
	    s->neighbor_accum_weight[eid] = sum;
	}
	s->sum_weight[nid] = sum;
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
  Version that can be used in synchronous or batch mode, where certain that node weights are already valid.
  And have already computed sum of weights for each node, and cumulative weight for each neighbor
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
#if DEBUG
    outmsg("Computing rat %d: node %d-->%d (%.3f/%.3f)", r, nid, g->neighbor[estart+offset], val, tsum);
#endif
    return g->neighbor[estart + offset];
}

/* Process single batch */
// TODO: Here's where things get interesting!
//    * Process rats currently in this zone
//    * Exchange rats (use function exchange_rats):
//      - Export rats that move out of this zone
//      - Import rats that move into this zone
//    * Exchange node states (use function exchange_node_states)
//      - Export counts for internal nodes adjacent to other zones
//      - Import counts for external nodes adjacent to this zone
//    * Compute weights for nodes in this zone
//    * Exchange weights (use function exchange_node_weights)
//      - Export weights for internal nodes adjacent to other zones
//      - Import weights for external nodes adjacent to this zone
static inline void do_batch(state_t *s, int batch, int bstart, int bcount) {
#if MPI
    int rank, size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank ); 
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    find_all_sums(s);
    START_ACTIVITY(ACTIVITY_NEXT);
    graph_t *g = s->g;
    int *zone_id = g->zone_id;
    int total_moveout=0;
    int *cur_len = s->cur_len;
    int *cur_rat = s->cur_rat;
    int *cur_rat_buff = s->cur_rat_buff;
    int buff_index = 0;
    int batch_size = s->batch_size;
    memset(s->moveout_count, 0, size*sizeof(int));

    int rid, onid, nnid;
    int batch_id = bstart/batch_size;
    for (int i=0;i<cur_len[batch_id];i++) {
        rid = cur_rat[batch_id*batch_size+i];
        onid = s->rat_position[rid];
        nnid = fast_next_random_move(s, rid);
        s->rat_position[rid] = nnid;
        s->rat_count[onid] -= 1;
        s->rat_count[nnid] += 1;

        if (zone_id[nnid]!=g->this_zone) {
            s->moveout_count[zone_id[nnid]]++;
            s->moveout_rat[total_moveout++] = rid;
        } else {
            cur_rat_buff[buff_index++] = rid;
        }
    }

    s->cur_len[batch_id] = buff_index;
    for (int i=0;i<buff_index;i++) {
        cur_rat[batch_id*batch_size+i] = cur_rat_buff[i];
    } 
    s->total_moveout = total_moveout;

    FINISH_ACTIVITY(ACTIVITY_NEXT);
    /* Update weights */
    exchange_rats(s);
    exchange_node_states(s);

    compute_cur_weights(s);
    exchange_node_weights(s);
#else
    int ri;
    find_all_sums(s);
    START_ACTIVITY(ACTIVITY_NEXT);
    for (ri = 0; ri < bcount; ri++) {
	int rid = ri+bstart;
	int onid = s->rat_position[rid];
	int nnid = fast_next_random_move(s, rid);
	s->rat_position[rid] = nnid;
	s->rat_count[onid] -= 1;
	s->rat_count[nnid] += 1;
    }
    FINISH_ACTIVITY(ACTIVITY_NEXT);
    /* Update weights */
    compute_all_weights(s);
#endif
}

static void batch_step(state_t *s) {
    int bstart = 0;
    int bsize = s->batch_size;
    int nrat = s->nrat;
    int bcount;
    int batch = 0;
    while (bstart < nrat) {
	bcount = nrat - bstart;
	if (bcount > bsize)
	    bcount = bsize;
	do_batch(s, batch, bstart, bcount);
	batch++;
	bstart += bcount;
    }
}

double simulate(state_t *s, int count, int dinterval, bool display) {
    int i;
    /* Compute and show initial state */
    bool show_counts = true;
    double start = currentSeconds();

#if MPI
    init_zone(s, s->g->this_zone);
#endif
    take_census(s);
    compute_all_weights(s);
    if (display) {
#if MPI
	if (s->g->this_zone == 0)
	    // Process 0 has a copy of the initial counts for all nodes.
	    show(s, show_counts);
#else
	show(s, show_counts);
#endif
    }
    for (i = 0; i < count; i++) {
	batch_step(s);
	if (display) {
	    show_counts = (((i+1) % dinterval) == 0) || (i == count-1);
#if MPI
	    if (s->g->this_zone == 0) {
		// Process 0 needs to call function show on each simulation step.
		// When show_counts is true, it will need to have
		// the counts for all other zones.
		// These must be gathered from the other processes.
		if (show_counts)
		    gather_node_state(s);
		show(s, show_counts);
	    } else {
		if (show_counts)
		    // Send counts to process 0
		    send_node_state(s);
	    }
#else
	    show(s, show_counts);
        
#endif
	}
    }
    double delta = currentSeconds() - start;
    done(s);
    return delta;
}

