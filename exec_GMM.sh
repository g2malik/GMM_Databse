#!/bin/bash
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=26
#SBATCH --time=01:00:00

module load gnu-parallel
module load python


# Run multiple cases at different snapshots

parallel -j 26 "python ./read_data_slices_SLURM.py 0 {}" ::: {52..1}