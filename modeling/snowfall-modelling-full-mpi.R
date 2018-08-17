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


dir_bm<-paste0(dir_R,"/00_biomod-full")
dir.create(dir_bm, showWarnings = FALSE, recursive = TRUE)
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
presmodstack_EA<-stack(paste0(dir_stacks,"present_modstack_EA.grd"))
r<-raster(nrow=nrow(presmodstack_EA),ncol=ncol(presmodstack_EA))
extent(r)<-extent(presmodstack_EA)
res(r)<-1000
crs(r)<-crs(presmodstack_EA)


presmodstack_EA<-resample(presmodstack_EA,r,method="bilinear")

#writeRaster( presmodstack_EA,paste0(dir_stacks, "/present_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
library(raster)

presmodstack<-stack(paste0(dir_stacks,"present_modstack_EA1km.grd"))


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
    setwd(dir_bm) 
    
    # # start log
    # get cat of consule and write to log file
    con <- file(paste(getwd(),"/",sp.n,"/",sp.n,".log", sep="") , open = "wt")
    
    # start sink log file
    sink(con, append=TRUE)
    sink(con, append=TRUE, type = "message")
   

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
    
    # retain only presence records (to generate psuedo-absence in BIOMOD_FormatingData)
    resp.occ.id <- which(pa[, myRespName] == 1)
    myResp <- as.numeric(pa[resp.occ.id, myRespName])
    
    # extract GuloGulo presences coordinates
    myRespXY <- pa[myResp,c('Longitude.x.','Latitude.y.')]
    
    # assign stacks to objects
    myExpl<-presmodstack

    # assign model metrics
    metrics = c( 'KAPPA', 'TSS', 'ROC')
    
    #################################################################
    # FORMAT INPUT DATA -- NOTE pseudo absences rep, selection, and #
    #################################################################
    
    # format input data for biomod
    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                         resp.xy = myRespXY,
                                         expl.var = myExpl,
                                         resp.name = myRespName,
                                         na.rm=TRUE,
                                         PA.nb.rep = 3,
                                         PA.nb.absences = 15000,
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
        NbRunEval=10,
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
    
    print(paste0("Building Ensemble Models for ",myRespName))
    
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
    
    print(paste0("Capturing Ensemble Model Outputs for  ",myRespName))
    
    # write em models built
    capture.output(get_built_models (myBiomodEM),
                   file=paste0(dir_out,"/",myRespName,"/",myRespName,"_em_models.txt"))
    
    # capture em model evals 
    print(paste0("Capturing Ensemble Models Evaluations ",myRespName))
    
    evalmodsEM<-get_evaluations(myBiomodEM,as.data.frame=TRUE)
    evalmodsEM$variety<-myRespName
    write.csv(evalmodsEM,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_models_EM_eval.csv"))
    
    #do.call(file.remove,list(list.files(pattern="temp*"))) 
    
    ### eval current model
    
    print(paste0("Capturing Model Ensemble Evaluations for ",myRespName))
    
    enevalmods<-get_evaluations(myBiomodEM,as.data.frame=TRUE)
    write.csv(enevalmods,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_em_evals-df.csv"))
    
    sink(type = "message")
    sink() 
    unlink(rtmpdir,recursive=TRUE)

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  

}
	


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

# # Stop the cluster
print("Stopping cluster")

sfStop()

#######
