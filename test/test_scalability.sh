#!/bin/bash

for p in 1 2; do
    for j in 1 2 4; do
        for k in 4 5 6 7 8; do
            n=$((2**k))
            echo "Running with $p MPI processes, $j openMP threads and grid size $n takes"
            mpirun -n $p ./main $n 4 test
        done
    done
done
