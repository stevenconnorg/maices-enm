source("init.R")

#install.packages("biomod2", repos = "http://r-forge.r-project.org",dependencies=TRUE)
library(biomod2)
library(raster)
library(rgdal)
library(biomod2)


#################################################################
# PRELIMINARY DATA FORMATTING
#################################################################



# read in model raster stacks in LONG/LAT
# can't project onto equal area because Maxent needs resolution with equal x, y
# this data came from lat/long, so the projected 

#presmodstack<-stack(paste0(dir_stacks,"present_modstack.grd"))
#f50modstack<-stack(paste0(dir_stacks,"f50_modstack.grd"))
#f70modstack<-stack(paste0(dir_stacks,"f70_modstack.grd"))

#equalarea <-
#  CRS(
#    "+proj=aea +lat_1=14.5 +lat_2=32.5 +lat_0=24 +lon_0=-105 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
#  )

#presmodstackEA <-
#  projectRaster(presmodstack, crs = equalarea, method = "bilinear")
#f50modstackEA <-
#  projectRaster(f50modstack, crs = equalarea, method = "bilinear")
#f70modstackEA <-
#  projectRaster(f70modstack, crs = equalarea, method = "bilinear")

#r<-raster(nrow=nrow(presmodstackEA),ncol=ncol(presmodstackEA))

# same extent
#extent(r)<-extent(presmodstackEA)

# but change to square resolution (1km^2)
#res(r)<-1000

# apply sample crs
#crs(r)<-crs(presmodstackEA)

## use bilinear if your predictors are continuous
#presmodstackEA1km<-resample(presmodstackEA,r,method="bilinear")
#f50modstackEA1km<-resample(f50modstackEA,r,method="bilinear")
#f70modstackEA1km<-resample(f70modstackEA,r,method="bilinear")

## function to define the intersect of rasters
## remove cells that dont have data in all layers
intersect_mask <- function(x){
  values_x <- raster::getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- raster::setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}


#presmodstackEA1km<- stack(mask(presmodstackEA1km, intersect_mask(presmodstackEA1km)))
#print("on future 50:")

#f50modstackEA1km<- stack(raster::mask(f50modstackEA1km, intersect_mask(f50modstackEA1km)))
#print("on future 70 stack:")

#f70modstackEA1km<- stack(raster::mask(f70modstackEA1km, intersect_mask(f70modstackEA1km)))


#writeRaster( presmodstackEA1km,paste0(dir_stacks, "/present_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
#writeRaster(f50modstackEA1km,paste0(dir_stacks, "/f50_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
#writeRaster(f70modstackEA1km,paste0(dir_stacks, "/f70_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)

presmodstackEA1km<-stack(paste0(dir_stacks,"present_modstack_EA1km.grd"))
f50modstackEA1km<-stack(paste0(dir_stacks,"f50_modstack_EA1km.grd"))
f70modstackEA1km<-stack(paste0(dir_stacks,"f70_modstack_EA1km.grd"))

presmodstack<-presmodstackEA1km
f50modstack<-f50modstackEA1km
f70modstack<-f70modstackEA1km


# check units
f50modstack[[1]]
f70modstack[[1]]
presmodstack[[1]]


#################################################################
# INITIALIZE FUNCTION TO APPLY TO EACH VARIETY
#################################################################

#allmodels<-c("GLM","GAM","GBM","ANN","CTA","RF","MARS","FDA","MAXENT.Phillips","MAXENT.Tsuruoka")
#models = c("MAXENT.Phillips")
# allmetrics = c(  'KAPPA', 'TSS', 'ROC', 'FAR','SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS')
metrics = c( 'KAPPA', 'TSS', 'ROC')


# following https://rpubs.com/dgeorges/190889

#maxent.background.dat.dir <- paste0(getwd(),"/maxent_bg")
#dir.create(maxent.background.dat.dir, showWarnings = FALSE, recursive = TRUE)


options(max.print=1000000)  # set max.print option high to capture outputs

setwd(dir_bmf) 

# get species vector name by getting all with ensemble projections for 2070 time period
# or whichever was last in the BioModApply function/loop

emprojfiles<-Sys.glob(paste0(dir_bmf,"/*/*/","*70.ensemble.projection.out")) # get ensemble files out
emproj.dirs<-sub(paste0(dir_bmf,"/"),"",emprojfiles) # remove everything from working directory path
sp.n<-sub(" */.*", "",emproj.dirs) # remove everything before first slash to get variety names that have ensemble models
sp.n


#for (sp.n in sp.n){
BioProjApply <-function(sp.n) {
  
  tryCatch({
    
    setwd(dir_bmf) 
    
    
    #con <- file(paste(getwd(),"/",sp.n,"/",sp.n,"proj.log", sep="") , open = "wt")
    #sink(con, append=TRUE)
    #sink(con, append=TRUE, type = "message")
    
    myRespName = sp.n
    
    myExpl<-presmodstack
    myExplFuture50<-f50modstack
    myExplFuture70<-f70modstack
    metrics = c( 'KAPPA', 'TSS', 'ROC')
    
    #load("/gpfs/home/scg67/thesis/02_R/maices-enm/.RData")
    #rasterOptions()$tmpdir      # get raster temp director
    #rtmpdir<-paste0(getwd(),"/tmp/",sp.n)
    #dir.create(rtmpdir, showWarnings = FALSE, recursive = TRUE)
    #rasterOptions(tmpdir=rtmpdir)  # set raster temp directory
    #rasterOptions(chunksize = 1e+07,maxmemory = 1e+09)
    
    
    myBiomodModelOut  <- get(load(paste0(getwd(),"/",myRespName,"/",myRespName,".",myRespName,"_current.models.out")))
    # 
    # myBiomodProj<-get(load(paste0(dir_bmf,"/",sp.n,"/proj_current/",sp.n,".current.projection.out")))

    library(rgdal)
    library(raster)
    
    
    # model projections
    myBiomodProj <- BIOMOD_Projection(
     modeling.output = myBiomodModelOut  ,
     new.env = myExpl,
     proj.name = 'current',
     selected.models = 'all',
     binary.meth = metrics,
     compress = TRUE,
     build.clamping.mask = FALSE
    )
    
    #print(paste0("Projecting  onto Future (2050) for ",myRespName))
    # future projections for rcp 85 period 50
    myBiomodProjFuture50 <- BIOMOD_Projection(
      modeling.output = myBiomodModelOut  ,
      new.env = myExplFuture50,
      proj.name = 'rcp85_50',
      selected.models = 'all',
      binary.meth = metrics,
      compress = TRUE,
      build.clamping.mask = FALSE
    )

    #print(paste0("Projecting onto Future (2070) for ",myRespName))
    # future projections for rcp 85 period 70
    myBiomodProjFuture70 <- BIOMOD_Projection(
      modeling.output = myBiomodModelOut  ,
      new.env = myExplFuture70,
      proj.name = 'rcp85_70',
      selected.models = 'all',
      binary.meth = metrics,
      compress = TRUE,
      build.clamping.mask = FALSE
    )

    #do.call(file.remove,list(list.files(pattern="temp*"))) 
    
    #myBiomodProj<-get(load(paste0(dir_bmf,"/",sp.n,"/proj_current/",sp.n,".current.projection.out")))    
    #myBiomodProjFuture50<-get(load(paste0(dir_bmf,"/",sp.n,"/proj_rcp85_50/",sp.n,".rcp85_50.projection.out")))
    #myBiomodProjFuture70<-get(load(paste0(dir_bmf,"/",sp.n,"/proj_rcp85_70/",sp.n,".rcp85_70.projection.out")))


    #################################################################
    # FORECAST EMSEMBLE MODELS BY CHOSEN METRICS
    #################################################################
    setwd(dir_bmf) 
    myBiomodEM <- get(load(paste0(getwd(),"/",myRespName,"/",myRespName,".",myRespName,"_currentensemble.models.out")))

    #print(paste0("Performing Ensemble Forcasting onto Current Data for ",myRespName))
    
    # current ensemble projection
    myBiomodEF <- BIOMOD_EnsembleForecasting(
      EM.output = myBiomodEM ,
      projection.output = myBiomodProj,
      selected.models = 'all',
     em.by='all',
      binary.meth=metrics,
      compress=TRUE)
    
    #print(paste0("Performing Ensemble Forcasting onto Future (2050) Data for ",myRespName))
    
    f50BiomodEF <- BIOMOD_EnsembleForecasting(
      EM.output = myBiomodEM ,
      projection.output = myBiomodProjFuture50,
      selected.models = 'all',
      em.by='all',
      binary.meth=metrics,
      compress=TRUE)
    
    #do.call(file.remove,list(list.files(pattern="temp*"))) 
    
    #print(paste0("Performing Ensemble Forcasting onto Future (2070) Data for ",myRespName))
    
    f70BiomodEF <- BIOMOD_EnsembleForecasting(
      EM.output = myBiomodEM ,
      projection.output = myBiomodProjFuture70,
      selected.models = 'all',
      em.by='all',
      binary.meth=metrics,
      compress=TRUE)
    
    #cat("\n\nExporting Ensemble as grd ...\n\n")
    
    #do.call(file.remove,list(list.files(pattern="temp*"))) 
    
    # alternatively, doing projections then ensemble by time-period/space
    ## do projections
    
    #cat("\n\nExporting Ensemble as grd ...\n\n")
    
    #do.call(file.remove,list(list.files(pattern="temp*"))) 
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

#save.image()


#######
# snowfall initialization
#install.packages("Rmpi", lib="/gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/", repos='http://cran.us.r-project.org')
library(snowfall)
library(Rmpi)
library(future)
library(parallel)

#args = commandArgs(trailingOnly = TRUE);
#ncpus = args[1];
# ------------------------------------------------------------------------
# initialize parallel mode using sockets and command-line args
# ------------------------------------------------------------------------
ntasks<-as.numeric(Sys.getenv("SLURM_NTASKS"))

print(ntasks)
print(class(ntasks))
print(paste0(ntasks," tasks"))

nnodes<-as.numeric(Sys.getenv("SLURM_JOB_NUM_NODES"))
print(nnodes)
print(class(nnodes))
print(paste0(nnodes," Nodes"))

ncpus<-ntasks*nnodes
#ncpus<-detectCores()
print(paste0(ncpus," CPUs"))

cl <- makeCluster(ncpus)
clusterEvalQ(cl, library(biomod2))
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(rgdal))


clusterExport(cl=cl, list("sp.n","BioProjApply","presmodstack","f50modstack","f70modstack",
                          "dir_out","dir_R","dir_dat","dir_stacks","dir_bmf","root"))

#for (sp.n in sp.n){
#print(sp.n)
parLapply(cl, sp.n[1:6],BioProjApply)
parLapply(cl, sp.n[7:15],BioProjApply)
parLapply(cl, sp.n[16:21],BioProjApply)
parLapply(cl, sp.n[22:29],BioProjApply)
parLapply(cl, sp.n[30:37],BioProjApply)
parLapply(cl, sp.n[38:46],BioProjApply)

#}
stopCluster(cl)#######