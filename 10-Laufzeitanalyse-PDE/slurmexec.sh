#!/bin/sh
> ids

for f in slurm/*.job
do
	sbatch $f >> ids
done
