#!/bin/sh
srun -N 3 -n 3 -p west -o output/12x1 mpirun -np 12 ./partdiff-par 1 2 512 2 2 512
srun -N 3 -n 3 -p west -o output/24x1 mpirun -np 24 ./partdiff-par 1 2 512 2 2 512
srun -N 3 -n 3 -p west -o output/1x12 mpirun -np 1 ./partdiff-par 12 2 512 2 2 512
srun -N 3 -n 3 -p west -o output/1x24 mpirun -np 1 ./partdiff-par 24 2 512 2 2 512
srun -N 3 -n 3 -p west -o output/2x6 mpirun -np 2 ./partdiff-par 6 2 512 2 2 512
srun -N 3 -n 3 -p west -o output/2x12 mpirun -np 2 ./partdiff-par 12 2 512 2 2 512
srun -N 3 -n 3 -p west -o output/12x2 mpirun -np 12 ./partdiff-par 2 2 512 2 2 512
