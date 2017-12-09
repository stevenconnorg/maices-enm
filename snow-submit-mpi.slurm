#!/bin/bash
#SBATCH --job-name=maices-enm
#SBATCH --nodes=4
#SBACH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --time=7-00:00:00
#SBATCH --partition=parallel
#SBATCH --mem=10000

module purge
module load R
module load openmpi/2.1.1

export LD_LIBRARY_PATH="/share/apps/openmpi/2.1.1/lib:$LD_LIBRARY_PATH"
export R_LIBS=/gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/lib
export R_LIBS=/share/apps/R/3.4.0/lib64/R/library

PBS_NODEFILE=`/usr/bin/generate_pbs_nodefile`

/bin/cat $PBS_NODEFILE

# How many cores total do we have?
#NO_OF_CORES=`cat $PBS_NODEFILE | egrep -v '^#'\|'^$' | wc -l | awk '{print $1}'`

# setup modules
#. /etc/profile.d/modules.sh

ulimit -c unlimited

# cp /gpfs/home/scg7/thesis/02_R/maxent -R $HOME

mpirun -np 1 R --vanilla --slave -f snow-modelling-slim.R 

#mpirun -x LD_LIBRARY_PATH -n 1 -hostfile $PBS_NODEFILE R CMD BATCH snow-modelling-slim.R

#R --slave -f snow-modelling.R?
#mpirun -x LD_LIBRARY_PATH -n 1 Rscript snow-modelling.R > snow-%j.out
