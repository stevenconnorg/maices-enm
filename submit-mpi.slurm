#!/bin/bash
#SBATCH --job-name=maices-enm
#SBATCH -p parallel
#SBATCH -N 65
#SBATCH --mail-type=end
#SBATCH --mail-user=scg67@txstate.edu
#SBATCH --exclusive
#SBATCH --time=7-00:00:00


module purge

module load R
module load openmpi/2.1.1

export LD_LIBRARY_PATH="/share/apps/openmpi/2.1.1/lib:$LD_LIBRARY_PATH"

export R_LIBS=/gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/lib
export R_LIBS=/share/apps/R/3.4.0/lib64/R/library

# We want a list of nodes, one per line.
# Use the PBS compatibility wrapper to make a PBS-style nodefile
PBS_NODEFILE=`/usr/bin/generate_pbs_nodefile`

/bin/cat $PBS_NODEFILE

# How many cores total do we have?
NO_OF_CORES=`cat $PBS_NODEFILE | egrep -v '^#'\|'^$' | wc -l | awk '{print $1}'`

# setup modules
. /etc/profile.d/modules.sh

ulimit -c unlimited



mpirun -n 1 -machinefile $PBS_NODEFILE Rscript  modelling.R	
