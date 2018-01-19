#!/bin/bash
#SBATCH --time=7-0:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=28
#SBATCH --mail-user=scg67@txstate.edu
#SBATCH --mail-type=END
#SBATCH --job-name=SNOWFALL
#SBATCH -p parallel
#SBATCH --exclusive

scontrol show job $SLURM_JOB_ID >> echo $SLURM_JOB_ID"scontrol_info.txt"

tic=`date +%s`
echo "Start Time = "`date`

cd $SLURM_SUBMIT_DIR
echo "working directory = " = pwd 


module purge
module load java/1.8.0_162
module load R
module load rocks-openmpi

# export JAVA_HOME="path that you found"
export $JAVA_HOME
export CLASSPATH=/gpfs/home/scg67/thesis/02_R/00_biomod/maxent.jar
export LD_LIBRARY_PATH=/opt/openmpi/lib:$LD_LIBRARY_PATH
export R_LIBS=/gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/lib
export R_LIBS=/share/apps/R/3.4.0/lib64/R/library

mkdir -p tmp
export TMPDIR=$(pwd)/tmp


echo "Launching snowfall script using R"

Rscript snowfall-modelling-mpi.R

echo "All Done!"

echo "End Time = "`date`
toc=`date +%s`

elapsedTime=`expr $toc - $tic`
echo "Elapsed Time = $elapsedTime seconds" 
