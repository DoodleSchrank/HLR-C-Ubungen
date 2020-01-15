#!/bin/sh

for f in slurm/*.job
do
	sbatch $f >> ids
done
