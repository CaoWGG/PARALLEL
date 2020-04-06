#include <stdlib.h>
#include <stdio.h>

#if MPI
#include <math.h>
#include <mpi.h>
#endif

#include "cycletimer.h"
#include "instrument.h"

/* Instrument different sections of program */
static char *activity_name[ACTIVITY_COUNT] = { "unknown", "startup", "compute_weights", "compute_sums", "find_moves", "local_comm", "global_comm"};

#if MPI
#define DATA_COUNT (ACTIVITY_COUNT+2)
#include "rutil.h"
#endif

static bool initialized = false;

static bool tracking = false;
static double global_start_time = 0.0;

#define MAXDEPTH 20

static activity_t activity_stack[MAXDEPTH];
static int stack_level = 0;

static double current_start_time = 0.0;

static double accum[ACTIVITY_COUNT];

void track_activity(bool enable) {
    tracking = enable;
}

static void init_instrument() {
    if (!tracking)
	return;
    if (initialized)
	return;
    initialized = true;
    global_start_time = currentSeconds();
    int a;
    for (a = 0; a < ACTIVITY_COUNT; a++) {
	accum[a] = 0.0;
    }
    stack_level = 0;
    activity_stack[stack_level] = ACTIVITY_NONE;
}

void start_activity(activity_t a) {
    if (!tracking)
	return;
    init_instrument();
    int olda = activity_stack[stack_level];
    double new_time = currentSeconds();
    accum[olda] += new_time - current_start_time;
    current_start_time = new_time;
    activity_stack[++stack_level] = a;
    if (stack_level >= MAXDEPTH) {
	fprintf(stderr, "Runaway instrumentation activity stack.  Disabling\n");
	tracking = false;
	return;
    }
}

void finish_activity(activity_t a) {
    if (!tracking)
	return;
    init_instrument();
    int olda = activity_stack[stack_level];
    if (a != olda) {
	fprintf(stderr, "Warning.  Started activity %s, but now finishing activity %s.  Disabling\n",
		activity_name[olda], activity_name[a]);
	tracking = false;
	return;
    }
    double new_time = currentSeconds();
    accum[olda] += (new_time - current_start_time);
    current_start_time = new_time;
    stack_level--;
    if (stack_level < 0) {
	fprintf(stderr, "Warning, popped off bottom of instrumentation activity stack.  Disabling\n");
	tracking = false;
	return;
    }
}

#if MPI
static void send_activity_data(int local_node_count, int local_edge_count) {
    double data[DATA_COUNT];
    int tag = 0;
    for (int a = 0; a < ACTIVITY_COUNT; a++)
	data[a] = accum[a];
    data[ACTIVITY_COUNT] = (double) local_node_count;
    data[ACTIVITY_COUNT+1] = (double) local_edge_count;
    MPI_Send(data, DATA_COUNT, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
}

static void get_activity_data(int zid, double *data) {
    int tag = 0;
    MPI_Recv(data, DATA_COUNT, MPI_DOUBLE, zid, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

#endif

void show_activity(FILE *f, int local_node_count, int local_edge_count) {
    if (!tracking)
	return;
    init_instrument();
    int a;
    double elapsed = currentSeconds() - global_start_time;
    double unknown = elapsed;
    for (a = 1; a < ACTIVITY_COUNT; a++)
	unknown -= accum[a];
    accum[0] = unknown;
#if MPI
    int nzone;
    int this_zone;
    int zid;
    MPI_Comm_size(MPI_COMM_WORLD, &nzone);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_zone);
    if (this_zone == 0) {
	double zdata[nzone][DATA_COUNT];
	double zelapsed[nzone];
	double rowdata[nzone];
	for (a = 0; a < ACTIVITY_COUNT; a++)
	    zdata[0][a] = accum[a];
	zdata[0][ACTIVITY_COUNT] = (double) local_node_count;
	zdata[0][ACTIVITY_COUNT+1] = (double) local_edge_count;
	zelapsed[0] = 0;
	for (zid = 1; zid < nzone ; zid++) {
	    get_activity_data(zid, zdata[zid]);
	    zelapsed[zid] = 0;
	}
	/* Print table of all data */
	fprintf(f, "Zone    ");
	for (zid = 0; zid < nzone; zid++)
	    fprintf(f, "%8d", zid);
	fprintf(f, "      Max      Mean  StdDev\n");
	fprintf(f, "Nodes   ");
	for (zid = 0; zid < nzone; zid++) {
	    rowdata[zid] = zdata[zid][ACTIVITY_COUNT];
	    fprintf(f, "%8d", (int) zdata[zid][ACTIVITY_COUNT]);
	}
	fprintf(f, " %8.0f  %8.1f%8.1f\n",
		data_max(rowdata, nzone), data_mean(rowdata, nzone), data_stddev(rowdata, nzone));
	fprintf(f, "Edges   ");
	for (zid = 0; zid < nzone; zid++) {
	    rowdata[zid] = zdata[zid][ACTIVITY_COUNT+1];
	    fprintf(f, "%8d", (int) zdata[zid][ACTIVITY_COUNT+1]);
	}
	fprintf(f, " %8.0f  %8.1f%8.1f\n",
		data_max(rowdata, nzone), data_mean(rowdata, nzone), data_stddev(rowdata, nzone));
	for (a = 0; a < ACTIVITY_COUNT; a++) {
	    fprintf(f, "        ");
	    for (zid = 0; zid < nzone; zid++) {
		zelapsed[zid] += zdata[zid][a];
		int ms = (int) (zdata[zid][a] * 1000.0);
		rowdata[zid] = ms;
		fprintf(f, "%8d", ms);
	    }
	    fprintf(f, " %8.0f  %8.1f%8.1f    %s\n",
		    data_max(rowdata, nzone), data_mean(rowdata, nzone), data_stddev(rowdata, nzone), activity_name[a]);
	}
	fprintf(f, "Elapsed ");
	for (zid = 0; zid < nzone; zid++) {
	    int ms = (int) (zelapsed[zid] * 1000.0);
	    rowdata[zid] = ms;
	    fprintf(f, "%8d", ms);
	}
	fprintf(f, " %8.0f  %8.1f%8.1f\n",
		data_max(rowdata, nzone), data_mean(rowdata, nzone), data_stddev(rowdata, nzone));
    } else {
	send_activity_data(local_node_count, local_edge_count);
    }
#else
    fprintf(f, "    %8d zones %8d edges\n", local_node_count, local_edge_count);
    for (a = 0; a < ACTIVITY_COUNT; a++) {
	if (accum[a] == 0.0)
	    continue;
	double ms = accum[a] * 1000.0;
	double pct = accum[a] / elapsed * 100.0;
	fprintf(f, "    %8d ms    %5.1f %%    %s\n", (int) ms, pct, activity_name[a]); 
    }
    fprintf(f, "    %8d ms    %5.1f %%    elapsed\n", (int) (elapsed * 1000.0), 100.0);
#endif
}
