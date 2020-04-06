#ifndef CRUN_H

/* Defining variable MPI enables use of MPI primitives */
#ifndef MPI
#define MPI 0
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#if MPI
#include <mpi.h>
#endif

/* Optionally enable debugging routines */
#ifndef DEBUG
#define DEBUG 0
#endif

/* Use statically defined ILD values */
#define STATIC_ILF 0

#if DEBUG
/* Setting TAG to some rat number makes the code track that rat's activity */
#define TAG 0
#endif

#include "rutil.h"
#include "cycletimer.h"
#include "instrument.h"


/*
  Definitions of all constant parameters.  This would be a good place
  to define any constants and options that you use to tune performance
*/

/* What is the maximum line length for reading files */
#define MAXLINE 1024

/* What is the batch size as a fraction of the number of rats */
#define BATCH_FRACTION 0.02

/* What is the base ILF */
#define BASE_ILF 1.75

/* What is the crossover between binary and linear search */
#define BINARY_THRESHOLD 4


/* Update modes */
typedef enum { UPDATE_SYNCHRONOUS, UPDATE_BATCH, UPDATE_RAT } update_t;

/* All information needed for graphrat simulation */

/* Parameter abbreviations
   N = number of nodes
   M = number of edge
   K = number of regions
   R = number of rats
   B = batch size
   Z = number of zones
 */

/* Representation of graph */
typedef struct {
    /* General parameters */
    int nnode;
    int nedge;
    int width;
    int height;
    int nzone;

    /* Graph structure representation */
    // Adjacency lists.  Includes self edge. Length=M+N.  Combined into single vector
    int *neighbor;
    // Starting index for each adjacency list.  Length=N+1
    int *neighbor_start;
    // For each node, zone identifier (number between 0 and Z-1).  Length=N
    int *zone_id;
#if STATIC_ILF
    // NOTE: This data removed.  ILFs are computed dynamically
    // Ideal load factor for each node.  (This value gets read from file but is not used.)  Length=N
    double *ilf;
#endif

    /**** Low-level details of a specific zone ****/
    int this_zone;
    /* How many nodes are in this zone */
    int local_node_count;
    /* How many edges are in this zone */
    int local_edge_count;
    /* Ordered list of nodes in this zone */
    int *local_node_list;
    /* For each other zone z, how many nodes in this zone have connections to nodes in z.  Length = Z */
    int *export_node_count;
    /* For each other zone z, lists of nodes in this zone with connections to nodes in z.  Length = Z */
    int **export_node_list;
    /* For each other zone z, how many nodes in z have connections to nodes in this zone.  Length = Z */
    int *import_node_count;
    /* For each other zone z, lists of nodes in z with connections to nodes in this zone.  Length = Z */
    int **import_node_list;

} graph_t;

/* Representation of simulation state */
typedef struct {
    graph_t *g;

    /* Number of rats */
    int nrat;


    /* Random seed controlling simulation */
    random_t global_seed;

    /* State representation */
    // Node Id for each rat.  Length=R
    int *rat_position;
    // Rat seeds.  Length = R
    random_t *rat_seed;

    /* Redundant encodings to speed computation */
    // Count of number of rats at each node.  Length = N.
    int *rat_count;
    // Store weights for each node.  Length = N
    double *node_weight;

    /* Computed parameters */
    double load_factor;  // nrat/nnnode
    int batch_size;      // Number of rats per batch

    // Memory to store sum of weights for each node's region.  Length = N
    double *sum_weight;
    // Memory to store cummulative weights for each node's region.  Length = M+N
    double *neighbor_accum_weight;

    // TODO: Declare the data structures you need to:
    // * Keep track of the rats in this zone
    // * Have storage for buffers you use to communicate with other zones
    int max_count;
    int *offset;
    int *buff;

    int *cur_rat;
    int *cur_rat_buff;
    int *cur_len;

    int total_moveout;
    int *moveout_count; // len=zone num
    int *moveout_rat; // zone num * nrat
    int *moveout_send;
    int *moveout_recv;

    int *rat_count_send;
    int *rat_count_recv;

    double *node_weight_send;
    double *node_weight_recv;

#if MPI
    MPI_Status status;
    MPI_Request *reqs;
#endif


} state_t;
    
/* Representation of a region.  Used by partitioner */
typedef struct {
    int id;
    int x;  // Left X
    int y;  // Upper Y
    int w;  // Width
    int h;  // Height
    int node_count;  // Number of nodes
    int edge_count;  // Number of (directed edges)
    int zone_id;     // Zone assigned by partitioner
} region_t;


/*** Function in partition.c ***/
/*
  This function should assign a zone id to every region in the graph.
  Zone IDs should range between 0 and Z-1.
  The updates should be made to the zone_id field of
  each entry in the region array.

  Note that each region entry contains information about the position,
  node count, and edge count of the region.
*/
void assign_zones(region_t *region_list, int nregion, int nzone);

/*** Functions in graph.c. ***/
graph_t *new_graph(int width, int height, int nedge, int nzone);

void free_graph();

graph_t *read_graph(FILE *gfile, int nzone);

#if DEBUG
void show_graph(graph_t *g);
#endif

#if MPI
void send_graph(graph_t *g);
graph_t *get_graph();
#endif

bool setup_zone(graph_t *g, int this_zone, bool verbose);
void clear_zone(graph_t *g);

/*** Functions in simutil.c ***/
/* Print message on stderr */
void outmsg(char *fmt, ...);

/* Allocate and zero arrays of int/double */
int *int_alloc(size_t n);
double *double_alloc(size_t n);

/* Read rat file and initialize simulation state */
state_t *read_rats(graph_t *g, FILE *infile, random_t global_seed);

/* Comparison function for qsort */
int comp_int(const void *ap, const void *bp);

/* Generate done message from simulator */
void done(state_t *s);

/* Print state of simulation */
/* show_counts indicates whether to include counts of rats for each node */
void show(state_t *s, bool show_counts);

/*** Functions in sim.c ***/

/* Run simulation.  Return elapsed time in seconds */
double simulate(state_t *s, int count, int dinterval, bool display);

/* Called after complete graph and all rats provided by master */
/* Return false if cannot allocate all required space */
bool init_zone(state_t *s, int zid);


#if MPI
/* These functions will support the required communication operations */

/* Called by process 0 to collect node states from all other processes */
void gather_node_state(state_t *s);

/* Called by other processes to send their node states to process 0 */
void send_node_state(state_t *s);

/* Called by process 0 to distribute rat state to all nodes */
void send_rats(state_t *s);

/* Called by other nodes to get rat state from master and set up state data structure */
state_t *get_rats(graph_t *g, random_t global_seed);

/* Move rats between zones as they migrate */
void exchange_rats(state_t *s);

/* Exchange node counts for boundary nodes between zones */
void exchange_node_states(state_t *s);

/* Exchange weights of nodes on boundaries */
void exchange_node_weights(state_t *s);

#endif // MPI


#define CRUN_H
#endif /* CRUN_H */
