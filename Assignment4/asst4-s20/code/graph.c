#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "crun.h"

graph_t *new_graph(int width, int height, int nedge, int nzone) {
    bool ok = true;
    graph_t *g = malloc(sizeof(graph_t));
    if (g == NULL)
	return NULL;
    int nnode = width * height;
    g->width = width;
    g->height = height;
    g->nnode = nnode;
    g->local_node_count = nnode;
    g->nedge = nedge;
    g->local_edge_count = nedge;
    g->nzone = nzone;
    g->neighbor = calloc(nnode + nedge, sizeof(int));
    ok = ok && g->neighbor != NULL;
    g->neighbor_start = calloc(nnode + 1, sizeof(int));
    ok = ok && g->neighbor_start != NULL;
#if STATIC_ILF
    g->ilf = calloc(nnode, sizeof(double));
    ok = ok && g->ilf != NULL;
#endif
    if (nzone > 0) {
	g->zone_id = calloc(nnode, sizeof(int));
	ok = ok && g->zone_id != NULL;
    } else
	g->zone_id = NULL;
    if (!ok) {
	outmsg("Couldn't allocate graph data structures");
	return NULL;
    }
    return g;
}

void free_graph(graph_t *g) {
    free(g->neighbor);
    free(g->neighbor_start);
#if STATIC_ILF
    free(g->ilf);
#endif
    free(g->zone_id);
    free(g);
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

/* Find ID of node given (x,y) coordinates */
static inline int find_node(graph_t *g, int x, int y) {
    // Row-major ordering
    return y * g->width + x;
}


/* Read in graph file and build graph data structure */
graph_t *read_graph(FILE *infile, int nzone) {
    char linebuf[MAXLINE];
    int width, height;
    int nnode, nedge;
    int nregion = 0;
    int i, hid, tid;
    double ilf;
    int nid, eid;
    int lineno = 0;

    // Read header information
    while (fgets(linebuf, MAXLINE, infile) != NULL) {
	lineno++;
	if (!is_comment(linebuf))
	    break;
    }
    if (sscanf(linebuf, "%d %d %d %d", &width, &height, &nedge, &nregion) < 3) {
	outmsg("ERROR. Malformed graph file header (line 1)\n");
	return NULL;
    }

    nnode = width * height;

    graph_t *g = new_graph(width, height, nedge, nzone);
    if (g == NULL)
	return g;

    nid = -1;
    // We're going to add self edges, so eid will keep track of all edges.
    eid = 0;
    for (i = 0; i < nnode; i++) {
	while (fgets(linebuf, MAXLINE, infile) != NULL) {
	    lineno++;
	    if (!is_comment(linebuf))
		break;
	}
	if (sscanf(linebuf, "n %lf", &ilf) != 1) {
	    outmsg("Line #%d of graph file malformed.  Expecting node %d\n", lineno, i+1);
	}
#if STATIC_ILF                        
	g->ilf[i] = ilf;
#endif
    }
    for (i = 0; i < nedge; i++) {
	while (fgets(linebuf, MAXLINE, infile) != NULL) {
	    lineno++;
	    if (!is_comment(linebuf))
		break;
	}
	if (sscanf(linebuf, "e %d %d", &hid, &tid) != 2) {
	    outmsg("Line #%d of graph file malformed.  Expecting edge %d\n", lineno, i+1);
	    return false;
	}
	if (hid < 0 || hid >= nnode) {
	    outmsg("Invalid head index %d on line %d\n", hid, lineno);
	    return false;
	}
	if (tid < 0 || tid >= nnode) {
	    outmsg("Invalid tail index %d on line %d\n", tid, lineno);
	    return false;
	}
	if (hid < nid) {
	    outmsg("Head index %d on line %d out of order\n", hid, lineno);
	    return false;
	    
	}
	// Starting edges for new node(s)
	while (nid < hid) {
	    nid++;
	    g->neighbor_start[nid] = eid;
	    // Self edge
	    g->neighbor[eid++] = nid;
	}
	g->neighbor[eid++] = tid;
    }
    while (nid < nnode-1) {
	// Fill out any isolated nodes
	nid++;
	g->neighbor[eid++] = nid;
    }
    g->neighbor_start[nnode] = eid;

    if (nregion > 0) {
	region_t *region_list = calloc(nregion, sizeof(region_t));
	if (region_list == NULL) {
	    fprintf(stderr, "Couldn't allocate space for region list\n");
	    return NULL;
	}
	for (i = 0; i < nregion; i++) {
	    int x, y, w, h;
	    while (fgets(linebuf, MAXLINE, infile) != NULL) {
		lineno++;
		if (!is_comment(linebuf))
		    break;
	    }
	    if (sscanf(linebuf, "r %d %d %d %d", &x, &y, &w, &h) != 4) {
		outmsg("Line #%d of graph file malformed.  Expecting region %d\n", lineno, i+1);
		return false;
	    }
	    region_list[i].id = i;
	    region_list[i].x = x; region_list[i].y = y; region_list[i].w = w; region_list[i].h = h;
	    region_list[i].node_count = w * h;
	    region_list[i].zone_id = 0;
	    int edge_count = 0;
	    /* Assign region IDs to nodes and compute number of edges  */
	    int dx, dy;
	    for (dx = x; dx < x+w; dx++)
		for (dy = y; dy < y+h; dy++) {
		    int nid = find_node(g, dx, dy);
		    edge_count += g->neighbor_start[nid+1] - g->neighbor_start[nid];
		}
	    region_list[i].edge_count = edge_count;
	}
	/* Use partitioning function to assign zones to nodes */
	if (nzone > 0) {
	    assign_zones(region_list, nregion, nzone);
	    /* Now go back through regions and assign zones to nodes */
	    for (i = 0; i < nregion; i++) {
		int x, y, w, h;
		x = region_list[i].x; y = region_list[i].y;
		w = region_list[i].w; h = region_list[i].h;
		int zid = region_list[i].zone_id;
		if (zid < 0 || zid > nzone) {
		    outmsg("Invalid zone number %d assigned to region %d.", zid, i);
		    return NULL;
		}
		int dx, dy;
		for (dx = x; dx < x+w; dx++)
		    for (dy = y; dy < y+h; dy++) {
			int nid = find_node(g, dx, dy);
			g->zone_id[nid] = zid;
		    }
	    }
	}
	free(region_list);
	outmsg("Loaded graph with %d nodes, %d edges, and %d regions, partitioned into %d zones \n", nnode, nedge, nregion, nzone);
    } else {
	outmsg("Loaded graph with %d nodes, %d edges, and %d regions\n", nnode, nedge, nregion);
    }
    //#if DEBUG
//    show_graph(g);
//#endif
    return g;
}

#if DEBUG
void show_graph(graph_t *g) {
    int nid, eid;
    outmsg("Graph\n");
    for (nid = 0; nid < g->nnode; nid++) {
	outmsg("%d:", nid);
	for (eid = g->neighbor_start[nid]; eid < g->neighbor_start[nid+1]; eid++) {
	    outmsg(" %d", g->neighbor[eid]);
	}
	outmsg("\n");
    }
    
}
#endif

#if MPI
/** MPI routines **/
void send_graph(graph_t *g) {
    /* Send basic graph parameters */
    int width = g->width;
    int height = g->height;
    int nedge = g->nedge;
    int nzone = g->nzone;
    int nnode = width * height;
    int params[4] = {width, height, nedge, nzone};
    MPI_Bcast(params, 4, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(g->neighbor, nedge+nnode, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(g->neighbor_start, nnode+1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(g->zone_id, nnode, MPI_INT, 0, MPI_COMM_WORLD);
}

graph_t *get_graph() {
    int params[4];
    MPI_Bcast(params, 4, MPI_INT, 0, MPI_COMM_WORLD);
    int width = params[0];
    int height = params[1];
    int nedge = params[2];
    int nzone = params[3];
    int nnode = width * height;
    graph_t *g = new_graph(width, height, nedge, nzone);
    if (g == NULL)
	return g;
    MPI_Bcast(g->neighbor, nedge+nnode, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(g->neighbor_start, nnode+1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(g->zone_id, nnode, MPI_INT, 0, MPI_COMM_WORLD);
    return g;
}
#endif

/* For verbose mode */
// Help in checking list generation
static void format_list(int *list, int count, char *buf) {
    char *pos = buf;
    sprintf(buf, "[");
    pos += 1;
    int i;
    for (i = 0; i < count; i++) {
	if (i >= 10) {
	    sprintf(pos, " ... ]");
	    return;
	}
	if (i > 0) {
	    sprintf(pos, ", ");
	    pos += 2;
	}
	int val = list[i];
	int len = val <= 9 ? 1 : val <= 99 ? 2 : 3;
	sprintf(pos, "%d", val);
	pos += len;
    }
    sprintf(pos, "]");
}

// Used to clear out information from one zone before setting up another
void clear_zone(graph_t *g) {
    g->local_node_count = 0;
    g->local_edge_count = 0;
    free(g->local_node_list); g->local_node_list = NULL;
    g->export_node_count = 0;
    free(g->export_node_list); g->export_node_list = NULL;
    g->import_node_count = 0;
    free(g->import_node_list); g->import_node_list = NULL;    
}

/* Set up zone-specific data structures */
/* Return false if something goes wrong */
bool setup_zone(graph_t *g, int this_zone, bool verbose) {
    g->this_zone = this_zone;
    int nzone = g->nzone;
    int nnode = g->nnode;
    int local_node_count = 0;
    int local_edge_count = 0;
    int nid, zid, idx;
    /* Allocate maximum possible nodes.  Trim this down later */
    g->local_node_list = calloc(nnode, sizeof(int));
    if (g->local_node_list == NULL) {
	outmsg("Couldn't allocate space for local nodes");
	return false;
    }

    g->export_node_count = calloc(nzone, sizeof(int));
    g->export_node_list = calloc(nzone, sizeof(int*));
    g->import_node_count = calloc(nzone, sizeof(int));
    g->import_node_list = calloc(nzone, sizeof(int*));
    if (g->export_node_count == NULL ||
	g->export_node_list == NULL ||
	g->import_node_count == NULL ||
	g->import_node_list == NULL) {
	outmsg("Couldn't allocate space for export/import info");
	return false;
    }

    /* Space for storing sets */
    /* Bit-vector representation of set of nodes */
    unsigned char *viewed_nodes = calloc(nnode, sizeof(unsigned char));
    if (viewed_nodes == NULL) {
	outmsg("Couldn't allocate space for node list");
	return false;
    }
    /* For each zone z, will allocate bit vector representation of subset of nodes in our zone */
    unsigned char **reverse_viewed_nodes = calloc(nzone, sizeof(unsigned char *));
    if (reverse_viewed_nodes == NULL) {
	outmsg("Couldn't allocate space for reverse node list");
	return false;
    }

    /*
      Pass one.  Over all nodes in graph, determine which ones are in this zone,
      and determine the import counts from other zones
    */
    for (nid = 0; nid < nnode; nid++) {
	if (g->zone_id[nid] != this_zone)
	    continue;
	g->local_node_list[local_node_count++] = nid;
	local_edge_count += (g->neighbor_start[nid+1] - g->neighbor_start[nid]);
	for (int eid = g->neighbor_start[nid]+1; eid < g->neighbor_start[nid+1]; eid++) {
	    int onid = g->neighbor[eid];
	    int ozid = g->zone_id[onid];
	    if (ozid != this_zone && !viewed_nodes[onid]) {
		g->import_node_count[ozid]++;
		viewed_nodes[onid] = 1;
	    }
	}
    }
    /* Now allocate data structures */
    g->local_node_count = local_node_count;
    g->local_edge_count = local_edge_count;
    g->local_node_list = realloc(g->local_node_list, local_node_count * sizeof(int));
    for (zid = 0; zid < nzone; zid++) {
	if (g->import_node_count[zid] > 0) {
	    reverse_viewed_nodes[zid] = calloc(local_node_count, sizeof(unsigned char));
	    g->import_node_list[zid] = calloc(g->import_node_count[zid], sizeof(int));
	    /* Allocate maximum size and adjust later */
	    g->export_node_list[zid] = calloc(local_node_count, sizeof(int));
	    if (reverse_viewed_nodes[zid]  == NULL ||
		g->import_node_list[zid] == NULL || g->export_node_list[zid] == NULL) {
		outmsg("Couldn't allocate space for export/import lists");
		return false;
	    }
	}
    }
    /* Pass two.  Iterate over local nodes.  Create import lists and export counts and lists */
    memset(viewed_nodes, 0, nnode * sizeof(unsigned char));
    memset(g->import_node_count, 0, nzone * sizeof(int));
    for (idx = 0; idx < local_node_count; idx++) {
	nid = g->local_node_list[idx];
	for (int eid = g->neighbor_start[nid]+1; eid < g->neighbor_start[nid+1]; eid++) {
	    int onid = g->neighbor[eid];
	    int ozid = g->zone_id[onid];
	    if (ozid != this_zone && !viewed_nodes[onid]) {
		g->import_node_list[ozid][g->import_node_count[ozid]++] = onid;
		viewed_nodes[onid] = 1;
	    }
	    if (ozid != this_zone && !reverse_viewed_nodes[ozid][idx]) {
		g->export_node_list[ozid][g->export_node_count[ozid]++] = nid;
		reverse_viewed_nodes[ozid][idx] = 1;
	    }
	}
    }
    /* Final cleanup.  Trim lists and dispose of temporary bit vectors.  Sort the import lists  */
    free(viewed_nodes);
    for (zid = 0; zid < nzone; zid++) {
	free(reverse_viewed_nodes[zid]);
	qsort(g->import_node_list[zid], g->import_node_count[zid], sizeof(int), comp_int);
	g->export_node_list[zid] = realloc(g->export_node_list[zid], g->export_node_count[zid] * sizeof(int));
    }
    free(reverse_viewed_nodes);
    // Verbose mode.  Print all the zone info.
    if (verbose) {
	char out_buf[100];
	format_list(g->local_node_list, g->local_node_count, out_buf);
	outmsg("Zone %d has %d nodes: %s", this_zone, g->local_node_count, out_buf);
	outmsg("Zone %d has %d edges", this_zone, g->local_edge_count);

	for (zid = 0; zid < nzone; zid++) {

	    if (g->export_node_count[zid] > 0) {
		format_list(g->export_node_list[zid], g->export_node_count[zid], out_buf);
		outmsg("Zone %d has %d nodes connected to zone %d: %s", this_zone, g->export_node_count[zid], zid, out_buf);
	    }
	    if (g->import_node_count[zid] > 0) {
		format_list(g->import_node_list[zid], g->import_node_count[zid], out_buf);
		outmsg("Zone %d has %d nodes in zone %d connected to it %s", this_zone, g->import_node_count[zid], zid, out_buf);
	    }
	}
    }
    return true;
}
