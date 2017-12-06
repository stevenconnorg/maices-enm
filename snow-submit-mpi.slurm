#!/bin/bash
#SBATCH --job-name=maices-enm
#SBATCH --nodes=4
#SBATCH --exclusive
#SBATCH --time=7-00:00:00
#SBATCH --partition=parallel

module purge
module load R
module load openmpi/2.1.1

export LD_LIBRARY_PATH="/share/apps/openmpi/2.1.1/lib:$LD_LIBRARY_PATH"
export R_LIBS=/gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/lib
export R_LIBS=/share/apps/R/3.4.0/lib64/R/library

mpirun -n 1 Rscript snow-modelling.R
