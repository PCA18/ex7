#!/usr/bin/env bash
#SBATCH --partition=bench
#SBATCH --nodes 1
for threads in 1 2 4 8 16 32
do
	srun ./relax.exe 200 $threads
done

		
