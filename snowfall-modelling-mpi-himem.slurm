#!/bin/bash
#SBATCH --time=7-0:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --mail-user=scg67@txstate.edu
#SBATCH --mail-type=END
#SBATCH --mem=900000
#SBATCH --job-name=SNOWFALL
#SBATCH -p himem
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --uid=scg67

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
export JAVA_HOME=$JAVA_HOME
export CLASSPATH=$CLASSPATH:/gpfs/home/scg67/thesis/02_R/00_biomod/maxent.jar
export CLASSPATH=$CLASSPATH:/gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/lib/dismo/java/maxent.jar
mkdir /gpfs/home/scg67/.java/.userPrefs

chmod -R 755 /gpfs/home/scg67/.java

export JAVA_OPTS="-Djava.util.prefs.systemRoot=/gpfs/home/scg67"
export JAVA_OPTS="-Djava.util.prefs.userRoot=/gpfs/home/scg67/"

chmod +x /gpfs/home/scg67/thesis/02_R/00_biomod/maxent.jar
chmod +x /gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/dismo/java/maxent.jar

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/home/scg67/gdal/2.2.2/lib
export GDAL_DATA=$LD_LIBRARY_PATH:/gpfs/home/scg67/gdal/2.2.2/share
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/home/scg67/proj-4.9.2/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/home/scg67/proj-4.9.2/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/home/scg67/proj-4.9.2/share/proj
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/home/scg67/lib

export R_LIBS=/gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/lib
export R_LIBS=/share/apps/R/3.4.0/lib64/R/library

chown -R scg67 /gpfs/home/scg67
chmod -R u+rwX /gpfs/home/scg67

#mkdir -p tmp
#export TMPDIR=$(pwd)/tmp


echo "Launching snowfall script using R"

echo "System Environment: "
env

Rscript snowfall-modelling-mpi.R

echo "All Done!"

echo "End Time = "`date`
toc=`date +%s`

elapsedTime=`expr $toc - $tic`
echo "Elapsed Time = $elapsedTime seconds" 

#echo "Run with R script"
cat snowfall-modelling-mpi.R

echo "Run with slurm script" 

cat snowfall-modelling-mpi.slurm

