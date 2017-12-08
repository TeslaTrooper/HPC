#!/bin/bash
for y in 4 8
do
    for x in 1024 2048 4096
    do
        for i in {1..5}
        do
            echo "Gebietsgröße: $x, NumThreads=$y"
            time OMP_NUM_THREADS=$y ./gameoflife $x $x 30
            echo -e "\n"
        done
    done 
done