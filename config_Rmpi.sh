module load openmpi
module load R

wget wget https //mpi4py.googlecode.com/files/mpi4py-1.3.1.tar.gz
gunzip -c openmpi-3.0.0.tar.gz | tar xf -
cd openmpi-3.0.0
./configure --disable-dlopen --enable-shared --prefix=/~/
make all install

# http://www.stats.uwo.ca/faculty/yu/Rmpi/install.htm
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/openmpi-3.0.0/lib
export MPI_ROOT=~/openmpi-3.0.0


R CMD INSTALL --configure-args="--with-Rmpi-type='OPENMPI' " Rmpi_0.6-6.tar.gz

also maybe worked:
R CMD INSTALL Rmpi_0.6-6.tar.gz --no-test-load


check R library paths
.libPaths()
