#!/bin/bash
#SBATCH --time=14-0:00:00
#SBATCH --nodes=8
#SBATCH --ntasks=28
#SBATCH --mail-user=scg67@txstate.edu
#SBATCH --mail-type=END
#SBATCH --job-name=SNOWFALL
#SBATCH -p parallel
#SBATCH --exclusive
#SBATCH --uid scg67

# get start time
tic=`date +%s`
echo "Start Time = "`date`


cd $SLURM_SUBMIT_DIR
cd ..
cd 00_biomod/
echo "working directory = " 
pwd

# create local tmp directory 
mkdir -p tmp
export TMPDIR=$(pwd)/tmp


module purge
module load java/1.8.0_162
module load R
module load rocks-openmpi

# export java home
export JAVA_HOME=$JAVA_HOME

# export .jar classpaths
export CLASSPATH=$CLASSPATH:/gpfs/home/scg67/thesis/02_R/00_biomod/maxent.jar
export CLASSPATH=$CLASSPATH:/gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/lib/dismo/java/maxent.jar

# if you don't have local user prefs
#mkdir /gpfs/home/scg67/.java/.userPrefs


echo "modifying java options"
chmod -R 755 /gpfs/home/scg67/.java
export JAVA_OPTS="-Djava.util.prefs.systemRoot=/gpfs/home/scg67"
export JAVA_OPTS="-Djava.util.prefs.userRoot=/gpfs/home/scg67/"
export JAVA_OPTS="-Djava.io.tmpdir=$TMPDIR"

echo "making maxent jars executable (dismo for BIOMOD_Tuning function)"
chmod +x /gpfs/home/scg67/thesis/02_R/00_biomod/maxent.jar
chmod +x /gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/dismo/java/maxent.jar

echo "exporting GDAL and associated libraries"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/home/scg67/gdal/2.2.2/lib
export GDAL_DATA=$LD_LIBRARY_PATH:/gpfs/home/scg67/gdal/2.2.2/share
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/home/scg67/proj-4.9.2/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/home/scg67/proj-4.9.2/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/home/scg67/proj-4.9.2/share/proj
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/home/scg67/lib

echo "export local R libraries"
export R_LIBS=/gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/lib
export R_LIBS=/share/apps/R/3.4.0/lib64/R/library

echo "own and modify local directories"
chown -R scg67 /gpfs/home/scg67
chmod -R u+rwX /gpfs/home/scg67



# get environment
echo "Env:"
env

# launch r script
echo "Launching snowfall script using R"
Rscript /gpfs/home/scg67/thesis/02_R/maices-enm/snowfall-modelling-mpi.R

echo "All Done!"

echo "End Time = "`date`
toc=`date +%s`

elapsedTime=`expr $toc - $tic`
echo "Elapsed Time = $elapsedTime seconds" 

echo "Run with"
cat snowfall-modelling-mpi.R
cat snowfall-modelling-mpi.slurm

cd tmp
rm -R *
