CC=gcc

MPI=-DMPI
MPICC = mpicc


DEBUG=0
INSTRUMENT=1
CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG) -DTRACK=$(INSTRUMENT) -std=gnu99
LDFLAGS= -lm
DDIR = ./data

CFILES = crun.c graph.c simutil.c sim.c rutil.c cycletimer.c instrument.c partition.c
HFILES = crun.h rutil.h cycletimer.h instrument.h

all: crun-seq crun-mpi


crun-seq: $(CFILES) $(HFILES) 
	$(CC) $(CFLAGS) -o crun-seq $(CFILES) $(LDFLAGS)

crun-mpi: $(CFILES) $(HFILES)
	$(MPICC) $(CFLAGS) $(MPI) -o crun-mpi $(CFILES) $(LDFLAGS)

demo1: grun.py
	@echo "Running Python simulator with text visualization.  Synchronous mode."
	./grun.py -g data/g-012x012-hlbrtX.gph -r data/r-012x012-c5.rats -n 10 -u s -v a -p 0.5

demo2: grun.py
	@echo "Running Python simulator with heat-map visualization.  Synchronous mode"
	./grun.py -g data/g-012x012-hlbrtX.gph -r data/r-012x012-c5.rats -n 10 -u s -v a -v h

demo3: crun-seq grun.py
	@echo "Batch mode. Note the smooth transitions"
	./crun-seq -g data/g-032x032-hlbrtZ.gph -r data/r-032x032-c10.rats  -n 20 | ./grun.py -d -v h

demo4: crun-seq grun.py
	@echo "Running on a 160x160 uniform graph"
	./crun-seq -g data/g-160x160-hlbrtA.gph -r data/r-160x160-c40.rats  -n 50 | ./grun.py -d -v h

demo5: crun-seq grun.py
	@echo "Running on a 160x160 fractal graph"
	./crun-seq -g data/g-160x160-hlbrtB.gph -r data/r-160x160-c40.rats  -n 50 | ./grun.py -d -v h

demo6: crun-seq grun.py
	@echo "Running on a 160x160 fractal graph, with rats initially along diagonal"
	./crun-seq -g data/g-160x160-hlbrtB.gph -r data/r-160x160-d40.rats  -n 50 | ./grun.py -d -v h

demo7: crun-seq grun.py
	@echo "Running on a 160x160 fractal graph, with rats initially distributed randomly"
	./crun-seq -g data/g-160x160-hlbrtB.gph -r data/r-160x160-r40.rats  -n 50 | ./grun.py -d -v h

demo8: crun-seq
	@echo "Demonstrating instrumentation.  No visualization"
	./crun-seq -g data/g-160x160-hlbrtB.gph -r data/r-160x160-c40.rats  -n 50 -q -I

clean:
	rm -f *~ *.pyc
	rm -rf *.dSYM
	rm -rf regression-cache check
	rm -f crun crun-seq crun-mpi
