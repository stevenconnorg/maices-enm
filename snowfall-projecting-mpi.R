#SLURMdat<-load("hpc.RData")
# establish directories
root<-"/gpfs/home/scg67/thesis"
#root<-"E:/thesis"
#setwd(root)
getwd()
# new directories for biomod

dir_dat<-paste0(root,"/01_data")
dir_R<-paste0(root,"/02_R")
dir_out<-paste0(root,"/03_output")
dir_figs<-paste0(root,"/04_figs")
dir_lit<-paste0(root,"/05_lit")
dir_comp<-paste0(root,"/06_comp")
dir_presentations<-paste0(root,"/07_pres")

dir_maices<-paste0(dir_dat,"/maices")
dir_ind<-paste0(dir_dat,"/ind")


dir_bm<-paste0(dir_R,"/00_biomod")
dir_topo<-paste0(dir_dat,"/topo")

# 
dir_clim<-paste0(dir_dat,"/clim")
dir_pres<-paste0(dir_clim,"/present")
dir_fut<-paste0(dir_clim,"/future")

dir_p.mosaics<-paste0(dir_pres,"/2.0/")
dir_f.mosaics<-paste0(dir_fut,"/1.4/")

dir_stacks<-paste0(dir_dat,"/stacks/")

#install.packages("biomod2", repos = "http://r-forge.r-project.org",dependencies=TRUE)

library(biomod2)
library(mgcv)


#################################################################
# PRELIMINARY DATA FORMATTING
#################################################################

setwd(dir_bm) 
ensemble.files<-Sys.glob(paste0(getwd(),"/*/*","*_currentensemble.models.out")) # get ensemble files out
ensemble.dirs<-sub(paste0(getwd(),"/"),"",ensemble.files) # remove everything from working directory path
sp.n<-sub(" */.*", "",ensemble.dirs) # remove everything before first slash to get variety names that have ensemble models
 


# read in model raster stacks with EQUAL AREA PROJECTION
library(raster)

presmodstack<-stack(paste0(dir_stacks,"present_modstack_EA.grd"))
f50modstack<-stack(paste0(dir_stacks,"f50_modstack_EA.grd"))
f70modstack<-stack(paste0(dir_stacks,"f70_modstack_EA.grd"))

names(presmodstack)
names(f50modstack)
names(f70modstack)

# plot(f50modstack)
# plot(presmodstack)
# plot(f70modstack)



#################################################################
# INITIALIZE FUNCTION TO APPLY TO EACH VARIETY
#################################################################

allmodels<-c("GLM","GAM","GBM","ANN","CTA","RF","MARS","FDA","MAXENT.Phillips","MAXENT.Tsuruoka")
models = c("MAXENT.Phillips")
metrics = c(  'KAPPA', 'TSS', 'ROC', 'FAR','SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS')


# following https://rpubs.com/dgeorges/190889

#maxent.background.dat.dir <- paste0(getwd(),"/maxent_bg")
#dir.create(maxent.background.dat.dir, showWarnings = FALSE, recursive = TRUE)


options(max.print=1000000)  # set max.print option high to capture outputs

BioModApply <-function(sp.n) {
  tryCatch({
    
    myRespName = sp.n
    #myRespName = sp.n[1]
    myResp <- as.numeric(pa[,myRespName])
    myRespXY = pa[,c('Longitude.x.','Latitude.y.')]
    
    myExpl<-presmodstack
    myExplFuture50<-f50modstack
    myExplFuture70<-f70modstack
    
    setwd(dir_bm) 
    bm_out_file <- load(paste0(getwd(),"/",myRespName,"/",myRespName,".",myRespName,"_current.models.out"))
    myBiomodModelOut <- get(bm_out_file)
    
    # model projections
     myBiomodProj <- BIOMOD_Projection(
     modeling.output = myBiomodModelOut,
     new.env = myExpl,
     proj.name = 'current',
     selected.models = 'all',
     binary.meth = metrics,
     filtered.meth = metrics,
     compress = TRUE,
     clamping.mask = T,
     output.format = '.grd')
    
    #print(paste0("Projecting onto Future (2070) for ",myRespName))
    # future projections for rcp 85 period 70
     myBiomodProjFuture70 <- BIOMOD_Projection(
       modeling.output = myBiomodModelOut,
       new.env = myExplFuture70,
       proj.name = 'rcp85_70',
       selected.models = 'all',
       binary.meth = metrics,
       filtered.meth = metrics,
       compress = TRUE,
       clamping.mask = T,
       output.format = '.grd')
    
     #print(paste0("Projecting  onto Future (2050) for ",myRespName))
    # future projections for rcp 85 period 50
     myBiomodProjFuture50 <- BIOMOD_Projection(
       modeling.output = myBiomodModelOut,
       new.env = myExplFuture50,
       proj.name = 'rcp85_50',
       selected.models = 'all',
       binary.meth = metrics,
       filtered.meth = metrics,
        compress = TRUE,
       clamping.mask = T,
       output.format = '.grd')
    
    #do.call(file.remove,list(list.files(pattern="temp*"))) 
    
    #################################################################
    # FORECAST EMSEMBLE MODELS BY CHOSEN METRICS
    #################################################################
    setwd(dir_bm) 
    bm_em.out_file <- load(paste0(getwd(),"/",myRespName,"/",myRespName,".",myRespName,"_currentensemble.models.out"))
    myBiomodEM <- get(bm_em.out_file)
    
    #print(paste0("Performing Ensemble Forcasting onto Current Data for ",myRespName))
    
    # current ensemble projection
    myBiomodEF <- BIOMOD_EnsembleForecasting(
      EM.output = myBiomodEM,
      projection.output = myBiomodProj,
      selected.models = 'all',
      binary.meth=metrics,
      filtered.meth=metrics,
      compress=TRUE)
    
    #print(paste0("Performing Ensemble Forcasting onto Future (2070) Data for ",myRespName))
    
    f70BiomodEF <- BIOMOD_EnsembleForecasting(
      EM.output = myBiomodEM,
      projection.output = myBiomodProjFuture70,
      selected.models = 'all',
      binary.meth=metrics,
      filtered.meth=metrics,
      compress=TRUE)
    
    #cat("\n\nExporting Ensemble as grd ...\n\n")
    
    #do.call(file.remove,list(list.files(pattern="temp*"))) 
    
    #print(paste0("Performing Ensemble Forcasting onto Future (2050) Data for ",myRespName))
    
    f50BiomodEF <- BIOMOD_EnsembleForecasting(
      EM.output = myBiomodEM,
      projection.output = myBiomodProjFuture50,
      selected.models = 'all',
      binary.meth=metrics,
      filtered.meth=metrics,
      compress=TRUE)
    
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

#save.image()


#######
# snowfall initialization
#install.packages("Rmpi", lib="/gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/", repos='http://cran.us.r-project.org')
library(snowfall)
library(Rmpi)
library(future)
#args = commandArgs(trailingOnly = TRUE);
#ncpus = args[1];

# ------------------------------------------------------------------------
# initialize parallel mode using sockets and command-line args
# ------------------------------------------------------------------------
ntasks<-as.numeric(Sys.getenv("SLURM_NTASKS"))

print(ntasks)
print(class(ntasks))
print(paste0(ntasks," tasks"))

ncputask<-as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
print(ncputask)
print(class(ncputask))
print(paste0(ncputask," N CPUs per task"))

ncpus<-ncputask*ntasks
#ncpus<-detectCores()
print(paste0(ncpus," CPUs"))

#ncpus<-future::availableCores()
#ncpus<-mpi.universe.size()-1
sfInit( parallel=TRUE, cpus=ncpus, type="MPI")
# Export packages to snowfall
sfLibrary('biomod2', character.only=TRUE)
sfExportAll()

# you may also use sfExportAll() to export all your workspace variables
## Do the run
mySFModelsOut <- sfLapply(sp.n,BioModApply)
#save(mySFModelsOut)
# stop snowfall

sfStop()
#######