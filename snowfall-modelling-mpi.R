sessionInfo()

#SLURMdat<-load("hpc.RData")
# establish directories

print("Setting Root Directory")
root<-"/gpfs/home/scg67/thesis"
#root<-"E:/thesis"
setwd(root)
getwd()
# new directories for biomod


print("Establishing directories")
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

#install.packages("biomod2", repos = "http://r-forge.r-project.org") #,dependencies=TRUE)

library(biomod2)
#library(mgcv)


#################################################################
# PRELIMINARY DATA FORMATTING
#################################################################

print("Reading in Presence/Absence Matrix")

### read in PAM with projected coordinates
pa<-data.frame(read.csv(file=paste0(dir_out,"/pa_dataframe_EA.csv")))

colnames(pa) 
pa_sp.n<-pa[,4:ncol(pa)]

# filter species by obs count
print("Filtering varieties with > 14 observations:")
i <- colSums(pa_sp.n,na.rm=T) > 14 
pa_sp.n<-pa_sp.n[,i]
colnames(pa_sp.n)

# drop unwanted columns by name to get vector of species names
#colnames(pa)
#drops <- c("X","Longitude.x.","Latitude.y.","X_Lambert", "Y_Lambert")
#pa_sp.n<-pa_sp.n[ , !(names(pa_sp.n) %in% drops)]

colnames(pa_sp.n)


names<-paste0(colnames(pa_sp.n))
sp.n= dput(names) # keep only species name, remove lat/long/etc. 

print("The species to be modeled:")
print("")
#sp.n

sp.n <- sub("_", ".", sp.n)
sp.n <- sub("_", ".", sp.n)
sp.n <- sub("_", ".", sp.n)

#sp.ngrep<-grep("\\.",sp.n)
#sp.n<-sp.n[sp.ngrep]
sp.n

 colnames(pa) <- sub("_", ".", colnames(pa))
 colnames(pa) <- sub("_", ".", colnames(pa))
 colnames(pa) <- sub("_", ".", colnames(pa))
 colnames(pa)


# read in model raster stacks with EQUAL AREA PROJECTION
library(raster)

print("Reading in Equal-Area modeling raster stack:")
#presmodstack_EA<-stack(paste0(dir_stacks,"present_modstack_EA.grd"))
#r<-raster(nrow=nrow(presmodstack_EA),ncol=ncol(presmodstack_EA))
#extent(r)<-extent(presmodstack_EA)
#res(r)<-1000
#crs(r)<-crs(presmodstack_EA)

#f50modstack_EA<-stack(paste0(dir_stacks,"f50_modstack_EA.grd"))
#f70modstack_EA<-stack(paste0(dir_stacks,"f70_modstack_EA.grd"))


#presmodstack_EA<-resample(presmodstack_EA,r,method="bilinear")
#f50modstack_proj<-resample(f50modstack_EA,r,method="bilinear")
#f70modstack_proj<-resample(f70modstack_EA,r,method="bilinear")
#presmodstack_proj<-presmodstack_EA


#writeRaster( presmodstack_proj,paste0(dir_stacks, "/present_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
#writeRaster(f50modstack_proj,paste0(dir_stacks, "/f50_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
#writeRaster(f70modstack_proj,paste0(dir_stacks, "/f70_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
library(raster)

presmodstack<-stack(paste0(dir_stacks,"present_modstack_EA1km.grd"))
#presmodstack_proj<-stack(paste0(dir_stacks,"present_modstack.grd"))
f50modstack<-stack(paste0(dir_stacks,"f50_modstack.grd"))
f70modstack<-stack(paste0(dir_stacks,"f70_modstack.grd"))
print("Getting stack names:")

#names<-names(presmodstack)
names
print("Running intersect_mask function:")

## function to define the intersect of rasters
intersect_mask <- function(x){
values_x <- raster::getValues(x)
inter_x <- values_x %*% rep(1,nlayers(x))
mask <- raster::setValues(subset(x,1),values = (inter_x>0))
return(mask)
}

## keep only all cells that are defined for all layers
print("on current stack:")
presmodstack<- stack(mask(presmodstack, intersect_mask(presmodstack)))
print("on future 50:")
#presmodstack_proj<- stack(mask(presmodstack_proj, intersect_mask(presmodstack_proj)))

print("on future 50:")
f50modstack<- stack(raster::mask(f50modstack, intersect_mask(f50modstack)))

print("on future 70 stack:")
f70modstack<- stack(raster::mask(f70modstack, intersect_mask(f70modstack)))

#writeRaster( presmodstack_proj,paste0(dir_stacks, "/present_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
#writeRaster(f50modstack,paste0(dir_stacks, "/f50_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
#writeRaster(f70modstack,paste0(dir_stacks, "/f70_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
library(raster)

#print("Reassigning stack names:")

#names(presmodstack)<-names
#names(f50modstack)<-names
#names(f70modstack)<-names


print("Reading in estimated prevalence .csv:")

prev<-as.data.frame(read.csv(file.path(dir_maices,"tau.csv")))
colnames(prev)
prev$Raza_prima
prev$Raza_prima <- sub(" ", ".", prev$Raza_prima)
prev$Raza_prima <- sub(" ", ".", prev$Raza_prima)
prev$Raza_prima <- sub(" ", ".", prev$Raza_prima)
prev$Raza_prima



#################################################################
# INITIALIZE FUNCTION TO APPLY TO EACH VARIETY
#################################################################


print("Establishing Function BioModApply() :")

setwd(dir_bm) 
options(max.print=1000000)  # set max.print option high to capture outputs

BioModApply <-function(sp.n) {
  tryCatch({
    
    # # start log
    # get cat of consule and write to log file
    #con <- file(paste(getwd(),"/",sp.n,"/",sp.n,".log", sep="") , open = "wt")
    
    # start sink log file
    #sink(con, append=TRUE)
    #sink(con, append=TRUE, type = "message")
   

    # # set some raster tmp file options

    # make tmp dir inside each sp.n dir
    rtmpdir<-paste0(getwd(),"/tmp/",sp.n)
    dir.create(rtmpdir, showWarnings = FALSE, recursive = TRUE)
    rasterOptions(tmpdir=rtmpdir)  # set raster temp directory

    # get raster tmp dir
    print("raster tmp dir:") 
    rasterOptions()$tmpdir     

    # increase or decrease these as appropriate
    #rasterOptions(chunksize = 1e+07,maxmemory = 1e+09)

    # get variety/species name
    myRespName = sp.n
    
    # get species X, Y coordinates
    myRespXY = pa[,c('Longitude.x.','Latitude.y.')]
    
    # get vector of 1's, 0's, or NA for presences, absences, and pseudo-absences that correspond to Long/Lat 
    myResp <- as.numeric(pa[,myRespName])

    # assign stacks to objects
    myExpl<-presmodstack
    myExplFuture50<-f50modstack
    myExplFuture70<-f70modstack

    # assign model metrics
    metrics = c( 'KAPPA', 'TSS', 'ROC','POD')
    
    #################################################################
    # FORMAT INPUT DATA -- NOTE pseudo absences rep, selection, and #
    #################################################################
    
    # format input data for biomod
    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                         resp.xy = myRespXY,
                                         expl.var = myExpl,
                                         resp.name = myRespName,
                                         na.rm=TRUE,
                                         PA.nb.rep =1,
                                         PA.nb.absences = 10000,
                                         PA.strategy = 'random')
    
    #################################################################
    # DEFINE MODEL OPTIONS
    #################################################################    

    # edit default options accordingly
    BIOMOD_ModelOptions <- BIOMOD_ModelingOptions(MAXENT.Phillips = list(memory_allocated = NULL,
                                                                         #background_data_dir = maxent.background.dat.dir, # https://rpubs.com/dgeorges/190889
                                                                         #maximumbackground = 10000,
                                                    maximumiterations = 10000,
                                                    visible = FALSE,
                                                    linear = TRUE,
                                                    quadratic = TRUE,
                                                    product = FALSE,
                                                    threshold = FALSE,
                                                    hinge = TRUE,
                                                    #lq2lqptthreshold = 80,
                                                    #l2lqthreshold = 10,
                                                    #hingethreshold = 10,
                                                    #beta_threshold = -1,
                                                    #beta_categorical = -1,
                                                    #beta_lqp = -1,
                                                    beta_hinge = -1,
                                                    betamultiplier = 2.5,
                                                    defaultprevalence = 0.1))

 
    #################################################################
    # BUILD MODELS
    #################################################################
    
    myBiomodModelOut <-
      BIOMOD_Modeling(
        myBiomodData, 
        models = c("MAXENT.Phillips"), 
        models.options = BIOMOD_ModelOptions, 
        NbRunEval=7,
        DataSplit=70, 
        Prevalence = prev$tau[prev$Raza_prima == myRespName], # species tau
        VarImport=10,
        models.eval.meth = metrics,
        SaveObj = TRUE,
        rescal.all.models = FALSE,
        do.full.models = FALSE,
        modeling.id = paste0(myRespName,"_current"))
    
    #################################################################
    # CAPTURE DATA INPUT AND MODEL OUTPUTS
    #################################################################
    
    # write data used for modelling
    capture.output(get_formal_data(myBiomodModelOut),
                   file=paste0(dir_out,"/",myRespName,"/",myRespName,"_model_data.txt"))

    # get model evaluations by metrics
    evalmods<-get_evaluations(myBiomodModelOut,as.data.frame=TRUE)
    evalmods$variety<-myRespName
    write.csv(evalmods,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_models_eval.csv"))
    
    ### get variable importance
    modevalimport<-get_variables_importance(myBiomodModelOut,as.data.frame=TRUE)
    colnames(modevalimport)

    write.csv(modevalimport,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_var_imp.csv"))
    

    # #################################################################
    # BUILD ENSEMBLE MODELS
    #################################################################
    
    #print(paste0("Building Ensemble Models for ",myRespName))
    
    # ensemble modeling
    myBiomodEM <- BIOMOD_EnsembleModeling(
      modeling.output = myBiomodModelOut,
      chosen.models = 'all',
      em.by='algo',
      eval.metric = metrics,
      eval.metric.quality.threshold = NULL,
      prob.mean = F,
      prob.cv =F,
      prob.ci = F,
      prob.ci.alpha = 0.05,
      prob.median = F,
      committee.averaging = F,
      prob.mean.weight = T,
      prob.mean.weight.decay = "proportional",
      VarImport = 10)
    
    #################################################################
    # CAPTURE ENSEMBLE MODEL OUTPUTS
    #################################################################
    
    #print(paste0("Capturing Ensemble Model Outputs for  ",myRespName))
    
    # write em models built
    capture.output(get_built_models (myBiomodEM),
                   file=paste0(dir_out,"/",myRespName,"/",myRespName,"_em_models.txt"))
    
    # capture em model evals 
    #print(paste0("Capturing Ensemble Models Evaluations ",myRespName))
    evalmodsEM<-get_evaluations(myBiomodEM,as.data.frame=TRUE)
    evalmodsEM$variety<-myRespName
    write.csv(evalmodsEM,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_models_EM_eval.csv"))
    
    #do.call(file.remove,list(list.files(pattern="temp*"))) 
    
    ### eval current model
    
    #print(paste0("Capturing Model Ensemble Evaluations for ",myRespName))
    enevalmods<-get_evaluations(myBiomodEM,as.data.frame=TRUE)
    write.csv(enevalmods,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_em_evals-df.csv"))
    


  ## do projections
  #proj_scen <- c("myExpl_proj", "myExplFuture50", "myExplFuture70")
  
  #for(scen in proj_scen){

    #sink(con, append=TRUE)
    #sink(con, append=TRUE, type = "message")
    #setwd(dir_bm) 


    #cat("\n> projections of ", scen)
 
    bm_out_file <- load(paste0(dir_bm,"/",myRespName,"/",myRespName,".",myRespName,"_current.models.out"))
    myBiomodModelOut <- get(bm_out_file)

    ## single model projections
    sp_proj <- BIOMOD_Projection(  modeling.output = myBiomodModelOut,
                                   new.env = get(scen),
                                   proj.name = scen,
                                   selected.models = 'all',
                                   binary.meth = metrics,
                                   filtered.meth = NULL,
                                   compress = TRUE,
                                   build.clamping.mask = FALSE,
                                   do.stack = FALSE )
    bm_em.out_file <- load(paste0(dir_bm,"/",myRespName,"/",myRespName,".",myRespName,"_currentensemble.models.out"))
    myBiomodEM <- get(bm_em.out_file)

    ## ensemble model projections
    sp_ens_proj <- BIOMOD_EnsembleForecasting(EM.output =  myBiomodEM,
                                              projection.output = sp_proj,
                                              binary.meth = metrics,
                                             compress = TRUE,
                                              do.stack = FALSE)
    #    sink(type = "message")
    #    sink() 
      }
  
     return(paste(sp.n," modelling completed !", sep=""))
  
    #sink(type = "message")
    #sink() 
    unlink(rtmpdir,recursive=TRUE)

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  

}
	


#save.image()


#myModelsOutlapply <-lapply(sp.n,BioModApply)

## alternatively with lapply
# myLapply_SFModelsOut <- lapply( sp.n, BioModApply)

# sph.umich.edu/biostat/computing/cluster/examples/r.html
#library(Rmpi)

# following http://www.glennklockwood.com/data-intensive/r/lapply-parallelism.html#3-1-lapply-halfway-to-parallel

# cl <- makeCluster( mpi.universe.size(), type="MPI" )

########
#or try...  following https://rcc.uchicago.edu/docs/software/environments/R/index.html
#np <- 27
#cl <- makeMPIcluster(np)
#worker.init <- function(packages) {
#  for (p in packages) {
#  library(p, character.only=TRUE) }
#NULL }
#clusterCall(cl, worker.init, c('biomod2','raster'))

#clusterExport(cl,c('sp.n','SLURMdat','presmodstack','f50modstack','f70modstack'))
#myModelsOut<-parLapply(cl,sp.n,BioModApply)
#mySFModelsOutCL <- clusterApply(cl, sp.n,BioModApply)

#stopCluster(cl)
#mpi.exit()
#########





#########
#library(Rmpi)
# https://rcc.uchicago.edu/docs/software/environments/R/index.html
# Initialize SNOW using MPI communication. The first line will get the number of
# MPI processes the scheduler assigned to us. Everything else is standard SNOW


#clusterExport(cl, SLURMdat)
#mySFModelsOut<-parLapplyLB(cl,sp.n,BioModApply)
#stopCluster(cl)
#mpi.exit()



# https://help.rc.ufl.edu/doc/R_MPI_Example
#ns <- mpi.universe.size() - 1
#mpi.spawn.Rslaves(nslaves=ns)
#cluster <- makeMPIcluster(np)

# Tell all slaves to return a message identifying themselves
#mpi.bcast.cmd( id <- mpi.comm.rank() )
#mpi.bcast.cmd( ns <- mpi.comm.size() )
#mpi.bcast.cmd( host <- mpi.get.processor.name() )
#mpi.remote.exec(paste("I am",mpi.comm.rank(),"of",mpi.comm.size()))

#mpi.parLapply(sp.n,BioModApply)


#mpi.close.Rslaves()
#mpi.quit()

########

#mpi.spawn.Rslaves(); cl = getMPIcluster()



## JOB 22645 -- WORKING??
########
#library(Rmpi)
#library(snow)
# https://www.osc.edu/~kmanalo/r_parallel
#slaves <- as.numeric(Sys.getenv(c("PBS_NP")))-1
#slaves <- mpi.universe.size()

#cl <- makeCluster(slaves, type = "MPI")
# or following http://homepage.divms.uiowa.edu/~luke/R/cluster/cluster.html
#cl <- getMPIcluster()

#clusterExport(cl, SLURMdat)
#library(snowfall)
#sfLibrary('biomod2', character.only=TRUE)
#sfExportAll()



#tick <- proc.time()
#mySFModelsOut<-parLapplyLB(cl,sp.n,BioModApply)
#tock <- proc.time() - tick

#cat("\nsnow w/ Rmpi test times using", slaves, "MPI slaves: \n")

#stopCluster(cl)
#mpi.quit()

########






########
#cl <- makeCluster(64, type = "MPI")
#clusterExport(cl, SLURMdat)

#mySFModelsOut<-clusterApply(cl,sp.n,fun=BioModApply)
#mySFModelsOut<-parLapplyLB(cl,sp.n,BioModApply)
#mySFModelsOut<-parLapply(cl,sp.n,BioModApply)
#stopCluster(cl)
#Rmpi:mpi.finalize()
#######




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

nnodes<-as.numeric(Sys.getenv("SLURM_JOB_NUM_NODES"))
print(nnodes)
print(class(nnodes))
print(paste0(nnodes," Nodes"))

ncpus<-ntasks*nnodes
#ncpus<-detectCores()
print(paste0(ncpus," CPUs"))



#library(parallel)
#cluster <- makeCluster(ncpus,type="MPI")

#print("Exporting environment to cluster")
#clusterExport(cl=cluster, list("sp.n","BioModApply","presmodstack","f50modstack","f70modstack","pa","prev","dir_out","dir_R","dir_dat","dir_stacks","dir_bm","root"),
#envir=environment())

#print("Exporting biomod2 package to cluster")
#clusterEvalQ(cluster, library(biomod2))

#print("Running >> parLapply(cluster, sp.n,BioModApply)") 
#applyOut<-parLapply(cluster, sp.n[1:5],BioModApply)
#save(applyOut)

#print("Stopping cluster") 
#stopCluster(cluster)






### snowfall ###

#ncpus<-future::availableCores()
#ncpus<-mpi.universe.size()-1
sfInit( parallel=TRUE, cpus=ncpus, type="MPI",
       slaveOutfile = paste0(Sys.getenv("SLURM_JOB_ID"),"_sfInit.log"))

# # Export packages to snowfall
print("Exporting packages to cluster")
sfLibrary('biomod2', character.only = TRUE)
sfLibrary('raster', character.only = TRUE)
sfLibrary('rgdal', character.only = TRUE)

# # Export environment to snowfall
print("Exporting environment to cluster")
#sfExportAll()

sfExport("sp.n","BioModApply","presmodstack","f50modstack","f70modstack",
         "pa","prev","dir_out","dir_R","dir_dat","dir_stacks","dir_bm","root"#,"maxent.background.dat.dir"
	)

# # Do the run
# use sfClusterApplyLB when machines have different specs
#print("Running >> sfLapply(sp.n,BioModApply)") # sfLapply
mySFModelsOut <- sfLapply(sp.n, BioModApply)

#mySFModelsOut <- sfLapply(sp.n[2:28], BioModApply)
#mySFModelsOut <- sfLapply(sp.n[21:40], BioModApply)
#mySFModelsOut <- sfLapply(sp.n[41:length(sp.n)], BioModApply)

# # Stop the cluster
print("Stopping cluster")

sfStop()

#######
