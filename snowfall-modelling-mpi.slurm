#!/bin/bash
#SBATCH --time=7-0:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=scg67@txstate.edu
#SBATCH --mail-type=END
#SBATCH --job-name=SNOWFALL
#SBATCH -p parallel
#SBATCH --exclusive
#SBATCH --mem-per-cpu=4096
tic=`date +%s`
echo "Start Time = "`date`

cd $SLURM_SUBMIT_DIR
echo "working directory = "$SLURM_SUBMIT_DIR


module purge
module load R
module load rocks-openmpi
#module load openmpi-1.8-x86_64

export LD_LIBRARY_PATH=/opt/openmpi/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH="/share/apps/intel/impi/2018.0.128/intel64/lib:$LD_LIBRARY_PATH"
export R_LIBS=/gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/lib
export R_LIBS=/share/apps/R/3.4.0/lib64/R/library

mkdir -p tmp
export TMPDIR=$(pwd)/tmp


# count the number of processors
#np=`srun hostname -s | wc -l`

# generate nodelist
#nodelist=`srun hostname -s | sort | tr '\n' ' '`

echo "Launching snowfall script using R"
#mpirun -np 1 R --vanilla --slave -f snowfall-modelling-mpi.R 
#srun R CMD BATCH snowfall-modelling-mpi.R
#R CMD BATCH --no-save --no-restore "--args $np $nodelist" snowfall-modelling.R
Rscript snowfall-modelling-mpi.R
echo "All Done!"

echo "End Time = "`date`
toc=`date +%s`

elapsedTime=`expr $toc - $tic`
echo "Elapsed Time = $elapsedTime seconds" 
