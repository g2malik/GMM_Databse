#!/bin/bash
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=01:00:00

module load gnu-parallel
module load python


# Run multiple cases at different snapshots

parallel -j 40 "python ./read_data_SLURM.py {}" ::: {500..10..10}