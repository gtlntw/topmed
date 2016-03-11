#!/bin/bash
#SBATCH -J tophat                  # A single job name for the array
#SBATCH -n 1                       # Number of cores
#SBATCH -N 1                       # All cores on one machine
#SBATCH -p serial_requeue          # Partition
#SBATCH --mem 4000                 # Memory request (4Gb)
#SBATCH -t 0-2:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o tophat_%A_%a.out        # Standard output
#SBATCH -e tophat_%A_%a.err        # Standard error
#SBATCH --array=1-30

 echo "${SLURM_ARRAY_TASK_ID}" > "${SLURM_ARRAY_TASK_ID}".txt 