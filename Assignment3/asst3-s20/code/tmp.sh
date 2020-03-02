#!/bin/bash

for i in 1 2 3 4 5 6 7 8 9 10 11 12
do
./benchmark.py -t $i -r 1 -I -b D > $i.out	
done
