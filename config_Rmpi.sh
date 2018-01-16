cd ~/sources
wget https://cran.r-project.org/src/contrib/Rmpi_0.6-6.tar.gz 
gunzip -c Rmpi_0.6-6.tar.gz | tar xf -

module load R

R CMD INSTALL Rmpi_0.6-6.tar.gz --configure-args=--with-mpi=/opt/openmpi

R
library(Rmpi)
