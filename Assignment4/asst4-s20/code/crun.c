/* C implementation of graphrats simulator */

#include <string.h>
#include <getopt.h>

#include "crun.h"

static void full_exit(int code) {
    done(NULL);
#if MPI
    MPI_Finalize();
#endif
    exit(code);
}

static void usage(char *name) {
#if MPI
    char *use_string = "-g GFILE -r RFILE [-n STEPS] [-s SEED] [-q] [-i INT] [-I]";
#else // !MPI
    char *use_string = "-g GFILE -r RFILE [-n STEPS] [-s SEED] [-q] [-i INT] [-I] [-z ZONE]";
#endif
    outmsg("Usage: %s %s\n", name, use_string);
    outmsg("   -h        Print this message\n");
    outmsg("   -g GFILE  Graph file\n");
    outmsg("   -r RFILE  Initial rat position file\n");
    outmsg("   -n STEPS  Number of simulation steps\n");
    outmsg("   -s SEED   Initial RNG seed\n");
    outmsg("   -q        Operate in quiet mode.  Do not generate simulation results\n");
    outmsg("   -i INT    Display update interval\n");
    outmsg("   -I        Instrument simulation activities\n");
#if !MPI
    outmsg("   -z ZONE   Test partitioning into ZONE zones without running simulation");
#endif
    full_exit(0);
}

int main(int argc, char *argv[]) {
    FILE *gfile = NULL;
    FILE *rfile = NULL;
    int steps = 1;
    int dinterval = 1;
    random_t global_seed = DEFAULTSEED;
    double secs = 0.0;
    int c;
    graph_t *g = NULL;
    state_t *s = NULL;
    bool instrument = false;
    bool display = true;
    bool show_zones_only = false;
    int process_count = 1;
    int this_zone = 0;
    int nzone = 0;
#if MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_zone);
    nzone = process_count;

#endif
    bool mpi_master = this_zone == 0;
#if MPI
    char *optstring = "hg:r:R:n:s:i:qI";
#else
    char *optstring = "hg:r:R:n:s:i:qIz:";
#endif
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch(c) {
        case 'h':
            if (!mpi_master) break;
            usage(argv[0]);
            break;
        case 'g':
            if (!mpi_master) break;
            gfile = fopen(optarg, "r");
            if (gfile == NULL) {
                outmsg("Couldn't open graph file %s\n", optarg);
		full_exit(1);
            }
            break;
        case 'r':
            if (!mpi_master) break;
            rfile = fopen(optarg, "r");
            if (rfile == NULL) {
                outmsg("Couldn't open rat position file %s\n", optarg);
                full_exit(1);
            }
            break;
        case 'n':
            steps = atoi(optarg);
            break;
        case 's':
            global_seed = strtoul(optarg, NULL, 0);
            break;
        case 'q':
            display = false;
            break;
        case 'i':
            dinterval = atoi(optarg);
            break;
        case 'I':
            instrument = true;
            break;
#if !MPI
	case 'z':
	    nzone = atoi(optarg);
	    show_zones_only = true;
	    break;
#endif 	    
        default:
            if (!mpi_master) break;
            outmsg("Unknown option '%c'\n", c);
            usage(argv[0]);
        }
    }

    TRACK_ACTIVITY(instrument);
    START_ACTIVITY(ACTIVITY_STARTUP);

    if (mpi_master) {
      	if (gfile == NULL) {
	    outmsg("Need graph file\n");
	    usage(argv[0]);
	}
	if (rfile == NULL && !show_zones_only) {
	    outmsg("Need initial rat position file\n");
	    usage(argv[0]);
	}
	g = read_graph(gfile, nzone);
	if (g == NULL) {
	    full_exit(1);
	}

	/* This only happens in sequential mode */
	if (show_zones_only) {
	    for (int z = 0; z < nzone; z++) {
		outmsg("*********** Setting up zone %d **********", z);
		if (!setup_zone(g, z, true)) {
		    full_exit(1);
		}
		clear_zone(g);
	    }
	    full_exit(0);
	}

	s = read_rats(g, rfile, global_seed);
	if (s == NULL) {
	    full_exit(1);
	}

#if MPI
        /* Master distributes the graph to the other processors */
	send_graph(g);
	if (!setup_zone(g, this_zone, false))
	    full_exit(1);
        /* Master distributes rats to the other processors */
	send_rats(s);
	// * Distribute copy of rats to other zones
	if (!init_zone(s, this_zone)) {
	    outmsg("Couldn't allocate space for zone %d data structures.  Exiting", this_zone);
	    full_exit(0);
        }
#endif
    } else {
	/* The other nodes receive the graph from the master */

#if MPI
	g = get_graph();
	if (g == NULL) {
	    outmsg("No graph.  Exiting");
	    full_exit(0);
	}
	if (!setup_zone(g, this_zone, false)) {
	    outmsg("Couldn't set up zones.  Exiting");
	    full_exit(0);
	}
	/* The other nodes receive the rats from the master */
	s = get_rats(g, global_seed);
	if (s == NULL) {
	    outmsg("No rats.  Exiting");
	    full_exit(0);
	}
	if (!init_zone(s, this_zone)) {
	    outmsg("Couldn't allocate space for zone %d data structures.  Exiting", this_zone);
	    full_exit(0);
        }

#endif
    }

    FINISH_ACTIVITY(ACTIVITY_STARTUP);

    if (mpi_master)
	outmsg("Running with %d processes.\n", process_count);
    // if (mpi_master)
	// Right now, run sequential simulator on master node
	// TODO: All processes should run simulator on their zones
	secs = simulate(s, steps, dinterval, display);
    if (mpi_master)
	outmsg("%d steps, %d rats, %.3f seconds\n", steps, s->nrat, secs);

    SHOW_ACTIVITY(stderr, g->local_node_count, g->local_edge_count);
#if MPI
    MPI_Finalize();
#endif    
    return 0;
}

