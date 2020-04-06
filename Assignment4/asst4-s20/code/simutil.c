#include "crun.h"

void outmsg(char *fmt, ...) {
#if MPI
    int process_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
    if (process_id != 0)
	fprintf(stderr, "Process %.2d|", process_id);
#endif    
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    bool got_newline = fmt[strlen(fmt)-1] == '\n';
    if (!got_newline)
	fprintf(stderr, "\n");
}

/* Allocate n int's and zero them out.  Maybe you could use multiple threads ... */
int *int_alloc(size_t n) {
    return (int *) calloc(n, sizeof(int));
}

/* Allocate n doubles's and zero them out.  Maybe you could use multiple threads ... */
double *double_alloc(size_t n) {
    return (double *) calloc(n, sizeof(double));
}

/* Allocate n random number seeds and zero them out.  */
static random_t *rt_alloc(size_t n) {
    return (random_t *) calloc(n, sizeof(random_t));
}

/* Allocate simulation state */
static state_t *new_rats(graph_t *g, int nrat, random_t global_seed) {
    int nnode = g->nnode;

    state_t *s = malloc(sizeof(state_t));
    if (s == NULL) {
	outmsg("Couldn't allocate storage for state\n");
	return NULL;
    }

    s->g = g;
    s->nrat = nrat;
    s->global_seed = global_seed;
    s->load_factor = (double) nrat / nnode;

    /* Compute batch size as max(BATCH_FRACTION * R, sqrt(R)) */
    int rpct = (int) (BATCH_FRACTION * nrat);
    int sroot = (int) sqrt(nrat);
    if (rpct > sroot)
	s->batch_size = rpct;
    else
	s->batch_size = sroot;

    // Allocate data structures
    bool ok = true;
    s->rat_position = int_alloc(nrat);
    ok = ok && s->rat_position != NULL;
    s->rat_seed = rt_alloc(nrat);
    ok = ok && s->rat_seed != NULL;
    s->rat_count = int_alloc(nnode);
    ok = ok && s->rat_count != NULL;

    s->node_weight = double_alloc(nnode);
    ok = ok && s->node_weight != NULL;
    s->sum_weight = double_alloc(g->nnode);
    ok = ok && s->sum_weight != NULL;
    s->neighbor_accum_weight = double_alloc(g->nnode + g->nedge);
    ok = ok && s->neighbor_accum_weight != NULL;

    if (!ok) {
	outmsg("Couldn't allocate space for %d rats", nrat);
	return NULL;
    }
    return s;
}

/* Set seed values for the rats.  Maybe you could use multiple threads ... */
static void seed_rats(state_t *s) {
    random_t global_seed = s->global_seed;
    int nrat = s->nrat;
    int r;
    for (r = 0; r < nrat; r++) {
	random_t seeds[2];
	seeds[0] = global_seed;
	seeds[1] = r;
	reseed(&s->rat_seed[r], seeds, 2);
#if DEBUG
	if (r == TAG)
	    outmsg("Rat %d.  Setting seed to %u\n", r, (unsigned) s->rat_seed[r]);
#endif
    }
}

/* See whether line of text is a comment */
static inline bool is_comment(char *s) {
    int i;
    int n = strlen(s);
    for (i = 0; i < n; i++) {
	char c = s[i];
	if (!isspace(c))
	    return c == '#';
    }
    return false;
}

/* Read in rat file */
state_t *read_rats(graph_t *g, FILE *infile, random_t global_seed) {
    char linebuf[MAXLINE];
    int r, nnode, nid, nrat;

    // Read header information
    while (fgets(linebuf, MAXLINE, infile) != NULL) {
	if (!is_comment(linebuf))
	    break;
    }
    if (sscanf(linebuf, "%d %d", &nnode, &nrat) != 2) {
	outmsg("ERROR. Malformed rat file header (line 1)\n");
	return false;
    }
    if (nnode != g->nnode) {
	outmsg("Graph contains %d nodes, but rat file has %d\n", g->nnode, nnode);
	return NULL;
    }
    
    state_t *s = new_rats(g, nrat, global_seed);


    for (r = 0; r < nrat; r++) {
	while (fgets(linebuf, MAXLINE, infile) != NULL) {
	    if (!is_comment(linebuf))
		break;
	}
	if (sscanf(linebuf, "%d", &nid) != 1) {
	    outmsg("Error in rat file.  Line %d\n", r+2);
	    return false;
	}
	if (nid < 0 || nid >= nnode) {
	    outmsg("ERROR.  Line %d.  Invalid node number %d\n", r+2, nid);
	    return false;
	}
	s->rat_position[r] = nid;
    }
    fclose(infile);

    seed_rats(s);
    outmsg("Loaded %d rats\n", nrat);
#if DEBUG
    outmsg("Load factor = %f\n", s->load_factor);
#endif
    return s;
}

/* print state of nodes */
void show(state_t *s, bool show_counts) {
    int nid;
    graph_t *g = s->g;
    printf("STEP %d %d %d\n", g->width, g->height, s->nrat);
    if (show_counts) {
	    for (nid = 0; nid < g->nnode; nid++)
		printf("%d\n", s->rat_count[nid]);
    }
    printf("END\n");
}

/* Print final output */
void done(state_t *s) {
#if MPI
    if (s == NULL || s->g->this_zone != 0)
	return;
#endif
    printf("DONE\n");
}

//TODO: Write function to initialize zone
bool init_zone(state_t *s, int zid) {
#if MPI
    int rank, size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank ); 
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    
    graph_t *g = s->g;
    int nrat = s->nrat;
    int nnode = g->nnode;
    int batch_size = s->batch_size;
    int *zone_id = g->zone_id;
    int *rat_position = s->rat_position;

    if (rank==0) {
        int *node_count = malloc(size*sizeof(int));
        MPI_Gather(&g->local_node_count,1,MPI_INT,node_count,1,MPI_INT,0,MPI_COMM_WORLD);
        int max_count = -1;
        for (int i=0;i<size;i++) {
            max_count = max_count>node_count[i]?max_count:node_count[i];
        }
        s->max_count = max_count;
        free(node_count);
    } else {
        MPI_Gather(&g->local_node_count,1,MPI_INT,NULL,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    MPI_Bcast(&(s->max_count),1,MPI_INT,0,MPI_COMM_WORLD);

    s->offset = malloc(size*sizeof(int));
    s->buff = malloc(size*s->max_count*sizeof(int));
    
    s->cur_rat = malloc(nrat*sizeof(int));
    s->cur_rat_buff = malloc(batch_size*sizeof(int));
    int nbatch = (nrat+batch_size-1)/batch_size, batch_id, id;
    s->cur_len= calloc(nbatch,sizeof(int));
    for (int rid=0;rid<nrat;rid++) {
        if (zone_id[rat_position[rid]]==zid) {
            batch_id = rid/batch_size;
            id = s->cur_len[batch_id];
            s->cur_rat[batch_id*batch_size+id] = rid;
            s->cur_len[batch_id]++;
        }
    }

    s->moveout_count = malloc(size*sizeof(int));
    s->moveout_rat = malloc(batch_size*sizeof(int));
    s->moveout_send = malloc(batch_size*3*sizeof(int));
    s->moveout_recv = malloc(batch_size*3*sizeof(int));

    s->rat_count_send = malloc(nnode*sizeof(int));
    s->rat_count_recv = malloc(nnode*sizeof(int));
    
    s->node_weight_send = malloc(nnode*sizeof(double));
    s->node_weight_recv = malloc(nnode*sizeof(double));

    s->reqs = malloc(size*sizeof(MPI_Request));


    MPI_Barrier(MPI_COMM_WORLD);

#endif
    return true;
}

//TODO: Implement these communication-support functions
#if MPI
/* Called by process 0 to distribute rat state to all nodes */
void send_rats(state_t *s) {
    START_ACTIVITY(ACTIVITY_GLOBAL_COMM);
    /* Your code should go here */
    MPI_Bcast(&s->nrat,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(s->rat_position,s->nrat,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(s->rat_seed,s->nrat,MPI_UNSIGNED,0,MPI_COMM_WORLD);
    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);
}

/* Called by other nodes to get rat state from master and set up state data structure */
state_t *get_rats(graph_t *g, random_t global_seed) {
    int nrat = 0;
    START_ACTIVITY(ACTIVITY_GLOBAL_COMM);
    /*
      Your code should go here.
      It should receive information about all rats, including how many there are,
      from process 0.
     */
    MPI_Bcast(&nrat,1,MPI_INT,0,MPI_COMM_WORLD);
    state_t *s = new_rats(g, nrat, global_seed);
    MPI_Bcast(s->rat_position,s->nrat,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(s->rat_seed,s->nrat,MPI_UNSIGNED,0,MPI_COMM_WORLD);
    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);
    return s;
}

void build_offset(int *offset, int *count, int size) {
    offset[0] = 0;
    for (int i=1;i<size;i++) {
        offset[i] = offset[i-1]+count[i-1];
    }
}

/* Called by process 0 to collect node states from all other processes */
void gather_node_state(state_t *s) {
    START_ACTIVITY(ACTIVITY_GLOBAL_COMM);
    /* Your code should go here */
    int rank, size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank ); 
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    graph_t *g = s->g;
    int max_count = s->max_count;
    int nnode = g->nnode;
    int *zone_id = g->zone_id;
    int *rat_count = s->rat_count;
    int *buff = s->buff;
    int *offset = s->offset;

    MPI_Gather(buff,max_count,MPI_INT,buff,max_count,MPI_INT,0,MPI_COMM_WORLD);
    memset(offset,0,size*sizeof(int));
    for (int i=0;i<nnode;i++) {
        int zid = zone_id[i];
        if (zid!=rank) {
            rat_count[i] = buff[zid*max_count+offset[zid]];
            offset[zid]+=1; 
        }
    }
    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);
}

/* Called by other processes to send their node states to process 0 */
void send_node_state(state_t *s) {
    START_ACTIVITY(ACTIVITY_GLOBAL_COMM);
    /* Your code should go here */
    graph_t *g = s->g;    
    int max_count = s->max_count;
    int *rat_count = s->rat_count;
    int *buff = s->buff;
    int local_node_count = g->local_node_count;
    int *local_node_list = g->local_node_list;
    for (int i=0;i<local_node_count;i++) {
        buff[i] = rat_count[local_node_list[i]];
    }
    MPI_Gather(buff,max_count,MPI_INT,NULL,0,MPI_INT,0,MPI_COMM_WORLD);
    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);
}

/* Move rats between zones as they migrate */
void exchange_rats(state_t *s) {
    START_ACTIVITY(ACTIVITY_COMM);
    /* Your code should go here */
    int rank, size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank ); 
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    graph_t *g = s->g;
    int *zone_id = g->zone_id;
    int *rat_count = s->rat_count;
    int *rat_position = s->rat_position;
    random_t *rat_seed = s->rat_seed;

    int *offset = s->offset;
    int total_moveout = s->total_moveout;
    int batch_size=s->batch_size;
    int *moveout_count = s->moveout_count; 
    int *moveout_rat = s->moveout_rat;
    int *moveout_send = s->moveout_send;
    int *moveout_recv = s->moveout_recv;

    MPI_Status status = s->status;
    MPI_Request *reqs = s->reqs;
    int req_cnt=0;

    for (int i=0;i<size;i++) {
        moveout_count[i] *= 3;
    }
    build_offset(offset, moveout_count, size);

    int zid, nid, rid, index;
    for (int i=0;i<total_moveout;i++) {
        rid = moveout_rat[i];
        nid = rat_position[rid];
        zid = zone_id[nid];
        index = offset[zid];
        moveout_send[index] = rid;
        moveout_send[index+1] = nid;
        moveout_send[index+2] = (int) rat_seed[rid];
        offset[zid]+=3;
    }

    // send 
    int start = 0;
    for (int i=0;i<size;i++) {
        if (i!=rank) {
            MPI_Isend(&moveout_send[start],moveout_count[i],MPI_INT,i,1,MPI_COMM_WORLD,reqs);
            start += moveout_count[i];
        }
    }

    for (int i=0;i<size;i++) {
        if (i!=rank) {
            MPI_Probe(i, 1, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_INT, &moveout_count[i]);
        }
    }

    build_offset(offset,moveout_count,size);

    req_cnt=0;
    for (int i=0;i<size;i++) {
        if (i!=rank) {
            MPI_Irecv(&moveout_recv[offset[i]],moveout_count[i],MPI_INT,i,1,MPI_COMM_WORLD,reqs+req_cnt);
            req_cnt+=1;
        }
    }
    MPI_Waitall(req_cnt,reqs,MPI_STATUS_IGNORE);
    
    int batch_id, id;
    for (int i=0;i<size;i++) {
        if (i!=rank) {
            int cnt = moveout_count[i]/3;
            for (int tmpid=0;tmpid<cnt;tmpid++) {
                int index = offset[i]+tmpid*3;
                rid = moveout_recv[index];
                nid = moveout_recv[index+1];
                rat_position[rid]=nid;
                rat_count[nid]+=1;
                rat_seed[rid]=(random_t)moveout_recv[index+2];

                batch_id=rid/batch_size;
                id=s->cur_len[batch_id];
                s->cur_rat[batch_id*batch_size+id] = rid;
                s->cur_len[batch_id]++;
            }
        }
    }
    FINISH_ACTIVITY(ACTIVITY_COMM);
}

/* Exchange node counts for boundary nodes between zones */
void exchange_node_states(state_t *s) {
    START_ACTIVITY(ACTIVITY_COMM);
    /* Your code should go here */
    int rank, size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank ); 
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    graph_t *g = s->g;
    int *export_node_count = g->export_node_count;
    int **export_node_list = g->export_node_list;
    int *import_node_count = g->import_node_count;
    int **import_node_list = g->import_node_list;
    int *rat_count = s->rat_count;

    int *rat_count_send = s->rat_count_send;
    int *rat_count_recv = s->rat_count_recv;
    int *offset = s->offset;

    MPI_Request *reqs = s->reqs;
    int req_cnt=0;

    build_offset(offset, export_node_count, size);

    int nid;
    for (int i=0;i<size;i++) {
        if (i!=rank) {
            for (int j=0;j<export_node_count[i];j++) {
                nid = export_node_list[i][j];
                rat_count_send[offset[i]+j] = rat_count[nid];
            }
        }
    }

    for (int i=0;i<size;i++) {
        if (i!=rank) {
            MPI_Isend(&rat_count_send[offset[i]],export_node_count[i],MPI_INT,i,0,MPI_COMM_WORLD,reqs);
            
        }
    }

    build_offset(offset,import_node_count,size);

    req_cnt = 0;
    for (int i=0;i<size;i++) {
        if (i!=rank) {
            MPI_Irecv(&(rat_count_recv[offset[i]]),import_node_count[i],MPI_INT,i,0,MPI_COMM_WORLD,reqs+req_cnt);
            req_cnt+=1;
        }
    }
    MPI_Waitall(req_cnt, reqs, MPI_STATUS_IGNORE);

    for (int i=0;i<size;i++) {
        if (i!=rank) {
            for (int j=0;j<import_node_count[i];j++) {
                nid = import_node_list[i][j];
                rat_count[nid] = rat_count_recv[offset[i]+j];
            }
        }
    }

    FINISH_ACTIVITY(ACTIVITY_COMM);
}


/* Exchange weights of nodes on boundaries */
void exchange_node_weights(state_t *s) {
    START_ACTIVITY(ACTIVITY_COMM);
    /* Your code should go here */

    int rank, size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank ); 
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    graph_t *g = s->g;
    int *export_node_count = g->export_node_count;
    int **export_node_list = g->export_node_list;
    int *import_node_count = g->import_node_count;
    int **import_node_list = g->import_node_list;
    double *node_weight = s->node_weight;

    double *node_weight_send = s->node_weight_send;
    double *node_weight_recv = s->node_weight_recv;
    int *offset = s->offset;

    MPI_Request *reqs = s->reqs;
    int req_cnt=0;

    build_offset(offset, export_node_count, size);

    int nid;
    for (int i=0;i<size;i++) {
        if (i!=rank) {
            for (int j=0;j<export_node_count[i];j++) {
                nid = export_node_list[i][j];
                node_weight_send[offset[i]+j] = node_weight[nid];
            }
        }
    }

    for (int i=0;i<size;i++) {
        if (i!=rank) {
            MPI_Isend(&node_weight_send[offset[i]],export_node_count[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,reqs);
            
        }
    }    
    
    build_offset(offset,import_node_count,size);

    req_cnt = 0;
    for (int i=0;i<size;i++) {
        if (i!=rank) {
            MPI_Irecv(&(node_weight_recv[offset[i]]),import_node_count[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,reqs+req_cnt);
            req_cnt+=1;
        }
    }
    MPI_Waitall(req_cnt, reqs, MPI_STATUS_IGNORE);
    
    for (int i=0;i<size;i++) {
        if (i!=rank) {
            for (int j=0;j<import_node_count[i];j++) {
                nid = import_node_list[i][j];
                node_weight[nid] = node_weight_recv[offset[i]+j];
            }
        }
    }    
    
    FINISH_ACTIVITY(ACTIVITY_COMM);
}

#endif // MPI

/* Function suitable for sorting arrays of int's */
int comp_int(const void *ap, const void *bp) {
    int a = *(int *) ap;
    int b = *(int *) bp;
    int lt = a < b;
    int gt = a > b;
    return -lt + gt;
}

