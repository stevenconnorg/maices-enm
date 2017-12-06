SLURMdat<-load("hpc.RData")
# establish directories
root<-"/gpfs/home/scg67/thesis"
#root<-"E:/thesis"
setwd(root)

# new directories for biomod


# dir_bmz<-paste0(dir_R,"/02_biomodez")


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

## recursively create directories if not already there 


## read in functions
# Function to Install and Load R Packages
# borrowed from Pratik Patil, https://stackoverflow.com/questions/9341635/check-for-installed-packages-before-running-install-packages
Install_And_Load <- function(Required_Packages)
{
  Remaining_Packages <- Required_Packages[!(Required_Packages %in% installed.packages()[,"Package"])];
  
  if(length(Remaining_Packages)) 
  {
    install.packages(Remaining_Packages,lib="home/scg67/R/x86_64_pc-linux-gnu-library/");
  }
  for(package_name in Required_Packages)
  {
    library(package_name,lib.loc="/home/scg67/R/x86_64_pc-linux-gnu-library/",character.only=TRUE,quietly=TRUE);
  }
}
# install and load required packages
requiredPackages<-(c("foreign",
                     "maptools",
                     "dplyr",
                     "rgdal",
                     "biomod2",
                     "reader",
                     "caret",
                     "tidyr",
                     "raster",
                     "quickPlot",
                     "biomod2",
                     "mgcv",
                     "gbm",
                     "dismo"
                   ))

library(biomod2)
library(mgcv)


#################################################################
# PRELIMINARY DATA FORMATTING
#################################################################

# set wd inside biomod subdirectory
setwd(dir_bm)


# get observation data formatted
pa<-read.csv(file=paste0(dir_out,"/pa_dataframe.csv"))
pa<-data.frame(pa)

names<-paste0(colnames(pa))

sp.n= dput(names [c(4:length(names))] # keep only species name, remove lat/long/etc. 
           ) #vector of species name(s), excluding lat and long cols


# read in model raster stacks

presmodstack<-stack(paste0(dir_stacks,"present_modstack.grd"))
f50modstack<-stack(paste0(dir_stacks,"f50_modstack.grd"))
f70modstack<-stack(paste0(dir_stacks,"f70_modstack.grd"))
# plot(f50modstack)
# plot(presmodstack)
# plot(f70modstack)




#################################################################
# INITIALIZE FUNCTION TO APPLY TO EACH VARIETY
#################################################################

allmodels<-c("GLM","GAM","GBM","ANN","CTA","RF","MARS","FDA","MAXENT.Phillips","MAXENT.Tsuruoka")
models = c("GLM","GAM","GBM")
metrics = c(  'KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS')



BioModApply <-function(sp.n) {

  load("/gpfs/home/scg67/thesis/02_R/maices-enm/.RData")
  rasterOptions()$tmpdir      # get raster temp directory
  rasterOptions(tmpdir=root)  # set raster temp directory
  options(max.print=1000000)  # set max.print option high to capture outputs
  maxentjar<-paste0(dir_R,"/maxent/maxent.jar") # define maxent jar location
  
  
  setwd(dir_bm) 

  myRespName = sp.n
  myResp <- as.numeric(pa[,myRespName])
  myRespXY = pa[,c('Longitude.x.','Latitude.y.')]
  
  myExpl<-presmodstack
  myExplFuture50<-f50modstack
  myExplFuture70<-f70modstack
  
  
  # Barbet-Massin et al 2012:
  #   Overall, we recommend the use of a large number (e.g. 10 000) of pseudo-absences with equal
  # weighting for presences and absences when using regression techniques (e.g. generalised linear
  # model and generalised additive model); averaging several runs (e.g. 10) with fewer pseudo-absences
  # (e.g. 100) with equal weighting for presences and absences with multiple adaptive regression splines
  # and discriminant analyses; and using the same number of pseudo-absences as available presences
  # (averaging several runs if few pseudo-absences) for classification techniques such as boosted regression
  # trees, classification trees and random forest. In addition, we recommend the random selection
  # of pseudo-absences when using regression techniques and the random selection of geographically
  # and environmentally stratified pseudo-absences when using classification and machine-learning
  # techniques
  
  
  #################################################################
  # FORMAT INPUT DATA -- NOTE pseudo absences rep, selection, and #
  #################################################################
  
  # format input data for biomod
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       resp.xy = myRespXY,
                                       expl.var = myExpl,
                                       resp.name = myRespName,
                                       PA.nb.rep = 10,
                                      PA.nb.absences = 500,
                                       PA.strategy = "sre",
										                  PA.sre.quant=0.025,
                                       na.rm=TRUE
  )
  
  
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  #################################################################
  # DEFINE MODEL OPTIONS
  #################################################################
  
  
  # print default biomod options 
  default_ModelOptions <-BIOMOD_ModelingOptions()
  print(default_ModelOptions)
  
  library(gam)
  # edit default options accordingly
  BIOMOD_ModelOptions <- BIOMOD_ModelingOptions(GLM = list( type = 'quadratic',      
                                                            interaction.level = 0,
                                                            myFormula = NULL,
                                                            test = 'AIC',           # tk BIC? 
                                                            family = binomial(link = 'logit'),
                                                            mustart = 0.5,
                                                            control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE) ),
                                                
                                                
                                                GBM = list( distribution = 'bernoulli',
                                                            n.trees = 2500,
                                                            interaction.depth = 7,
                                                            n.minobsinnode = 5,
                                                            shrinkage = 0.001,
                                                            bag.fraction = 0.5,
                                                            train.fraction = 1,
                                                            cv.folds = 3,
                                                            keep.data = FALSE,
                                                            verbose = FALSE,
                                                            perf.method = 'cv'),
                                                
                                                GAM = list(algo = 'GAM_mgcv', type = 's_smoother', k = NULL, 
                                                           interaction.level = 0, 
                                                           myFormula = NULL, 
                                                           family = 'binomial', 
                                                           control = gam.control(epsilon = 1e-06, trace = FALSE, maxit = 100)),
                                                
                                                CTA = list( method = 'class',
                                                            parms = 'default',
                                                            cost = NULL,
                                                            control = list(xval = 5, minbucket = 5, minsplit = 5, cp = 0.01, maxdepth = 25) ),
                                                
                                                
                                                ANN = list( NbCV = 5,
                                                            size = NULL,
                                                            decay = NULL,
                                                            rang = 0.1,
                                                            maxit = 200),
                                                
                                                SRE = list( quant = 0.025),
                                                
                                                FDA = list( method = 'mars',
                                                            add_args = NULL),
                                                
                                                MARS = list( type = 'simple',
                                                             interaction.level = 0,
                                                             myFormula = NULL,
                                                             nk = NULL,
                                                             penalty = 2,
                                                             thresh = 0.001,
                                                             nprune = NULL,
                                                             pmethod = 'backward'),
                                                
                                                RF = list( do.classif = TRUE,
                                                           ntree = 500,
                                                           mtry = 'default',
                                                           nodesize = 5,
                                                           maxnodes = NULL),
                                                
                                                MAXENT.Phillips = list( path_to_maxent.jar = maxentjar,
                                                                        memory_allocated = 1020,
                                                                        background_data_dir = 'default',
                                                                        maximumbackground = 'default',
                                                                        maximumiterations = 200,
                                                                        visible = FALSE,
                                                                        linear = TRUE,
                                                                        quadratic = TRUE,
                                                                        product = TRUE,
                                                                        threshold = TRUE,
                                                                        hinge = TRUE,
                                                                        lq2lqptthreshold = 80,
                                                                        l2lqthreshold = 10,
                                                                        hingethreshold = 15,
                                                                        beta_threshold = -1,
                                                                        beta_categorical = -1,
                                                                        beta_lqp = -1,
                                                                        beta_hinge = -1,
                                                                        betamultiplier = 1,
                                                                        defaultprevalence = 0.5),
                                                
                                                MAXENT.Tsuruoka = list( l1_regularizer = 0,
                                                                        l2_regularizer = 0,
                                                                        use_sgd = FALSE,
                                                                        set_heldout = 0,
                                                                        verbose = FALSE))
  
  
  #################################################################
  # TUNE MODEL OPTIONS -- NOT WORKING RN
  #################################################################
  
  # install.packages("ENMeval")
  # library(caret)
  # library(ENMeval)
  
  # download new version of code from Frank Breiner (writer of BIOMOD_tuning), attached here: http://r-forge.wu.ac.at/forum/forum.php?max_rows=75&style=nested&offset=152&forum_id=995&group_id=302
  
  # source(paste0(dir_R,"/maices-enm/BIOMOD.tuning_v6.R"))
  # library(doParallel);cl<-makeCluster(8);registerDoParallel(cl) 
  # devtools::install_github('topepo/caret/pkg/caret')
  # library(caret)
  # BIOMOD_TunedOptions <- BIOMOD_tuning(myBiomodData,
  #                                env.ME = myExpl,
  #                                n.bg.ME = ncell(myExpl)
  #                                )
  # stopCluster(cl)
  # BIOMOD_ModelOptions<-Biomod.tuning$models.options
  

  # capture.output(BIOMOD_TunedOptions$models.options,file=paste0(dir_out,"/model-opts/",myRespName,"_tuned_opts.txt"))
  
  
  #################################################################
  # BUILD MODELS
  #################################################################
  
  # modeling
  myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData, 
    models = models, 
    models.options = BIOMOD_ModelOptions, 
    NbRunEval=10,
    DataSplit=50,
    VarImport=5,
    models.eval.meth = metrics,
    SaveObj = TRUE,
    rescal.all.models = FALSE,
    do.full.models = TRUE,
    modeling.id = paste0(myRespName,"_current"))
  
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  print(paste0("Done Running Models for ",sp.n))
  
  dir.create(paste0(dir_out,"/",myRespName))
  
  

  #################################################################
  # CAPTURE DATA INPUT AND MODEL OUTPUTS
  #################################################################
  
  # write data used for modelling
  capture.output(get_formal_data(myBiomodModelOut),
                 file=paste0(dir_out,"/",myRespName,"/",myRespName,"_model_data.txt"))
  
  
  print(paste0("Capturing Model Evaluations for ",myRespName))
  evalmods<-get_evaluations(myBiomodModelOut,as.data.frame=TRUE)
  write.csv(evalmods,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_models_eval.csv"))
  
  ### get variable importance
  modevalimport<-get_variables_importance(myBiomodModelOut,as.data.frame=TRUE)
  write.csv(modevalimport,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_var_imp.csv"))
  
  ### get model summaries
  #capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"_current/",myRespName,"_PA1_Full_ANN",sep="")))))
                # ,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_ANN_summary.txt"))
  
  # capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"_current/",myRespName,"_PA1_Full_CTA",sep="")))))
                 # ,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_CTA_summary.txt"))
  
  # fda1<-get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"_current/",myRespName,"_PA1_Full_FDA",sep=""))))
  # capture.output(as.table(fda1$confusion),file==paste0(dir_out,"/",myRespName,"/",myRespName,"_FDA-conf_summary.csv"))

   # capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"_current/",myRespName,"_PA1_Full_GAM",sep="")))))
   #              ,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_GAM_summary.txt"))
    
  #capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"_current/",myRespName,"_PA1_Full_GBM",sep="")))))
  #               ,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_GBM_summary.txt"))
  
  #capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"_current/",myRespName,"_PA1_Full_GLM",sep="")))))
  #               ,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_GLM_summary.txt"))
  
  #capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"_current/",myRespName,"_PA1_Full_MARS",sep="")))))
  #               ,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_MARS_summary.txt"))
  
  # copy()
  # summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_PA1_Full_MAXENT.Phillips",sep="")))))
  # maxent_t1<-get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_PA1_Full_MAXENT.Tsuruoka",sep=""))))
  
  
  #################################################################
  # PROJECT MODELS ONTO CURRENT AND FUTURE CONDITIONS
  #################################################################
  
  print(paste0("Projecting onto Current Dataset for ",myRespName))
  
  # model projections
  myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExpl,
    proj.name = 'current',
    selected.models = models,
    binary.meth = metrics,
	filtered.meth = metrics,
    compress = TRUE,
    clamping.mask = T,
    output.format = '.grd')
  
  print(paste0("Projecting onto Future (2070) for ",myRespName))
  
  # future projections for rcp 85 period 70
  myBiomodProjFuture70 <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExplFuture70,
    proj.name = 'rcp85_70',
    selected.models = models,
    binary.meth = metrics,
	filtered.meth = metrics,
    compress = TRUE,
    clamping.mask = T,
    output.format = '.grd')
  
  print(paste0("Projecting  onto Future (2050) for ",myRespName))
  
  # future projections for rcp 85 period 50
  myBiomodProjFuture50 <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExplFuture50,
    proj.name = 'rcp85_50',
    selected.models = models,
    binary.meth = metrics,
	filtered.meth = metrics,
    compress = TRUE,
    clamping.mask = T,
    output.format = '.grd')
  
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  #################################################################
  # BUILD ENSEMBLE MODELS
  #################################################################
  
  print(paste0("Building Ensemble Models for ",myRespName))
  
  # ensemble modeling
  myBiomodEM <- BIOMOD_EnsembleModeling(
    modeling.output = myBiomodModelOut,
    chosen.models = models,
    em.by="algo",
    eval.metric = metrics,
    eval.metric.quality.threshold = c(rep(0.75,length(metrics))),
    prob.mean = T,
    prob.cv = T,
    prob.ci = F,
    prob.median = F,
    committee.averaging = T,
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
  capture.output(getEMeval(myBiomodEM),
                 file=paste0(dir_out,"/",myRespName,"/",myRespName,"_em_mods_eval.txt"))
  
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  ### eval current model
  
  print(paste0("Capturing Model Ensemble Evaluations for ",sp.n))
  enevalmods<-get_evaluations(myBiomodEM,as.data.frame=TRUE)
  write.csv(enevalmods,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_em_evals-df.csv"))
  
  #################################################################
  # FORECAST EMSEMBLE MODELS BY CHOSEN METRICS
  #################################################################
  
  print(paste0("Performing Ensemble Forcasting onto Current Data for ",myRespName))
  
  # current ensemble projection
  myBiomodEF <- BIOMOD_EnsembleForecasting(
    EM.output = myBiomodEM,
    projection.output = myBiomodProj,
    selected.models = models,
	binary.meth=metrics,
	filtered.meth=metrics,
	compress=TRUE)
  

  print(paste0("Performing Ensemble Forcasting onto Future (2070) Data for ",myRespName))
  
  f70BiomodEF <- BIOMOD_EnsembleForecasting(
    EM.output = myBiomodEM,
    projection.output = myBiomodProjFuture70,
    selected.models = models,
	binary.meth=metrics,
	filtered.meth=metrics,
	compress=TRUE)
  
  cat("\n\nExporting Ensemble as grd ...\n\n")
  
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  
  print(paste0("Performing Ensemble Forcasting onto Future (2050) Data for ",myRespName))
  
  f50BiomodEF <- BIOMOD_EnsembleForecasting(
    EM.output = myBiomodEM,
    projection.output = myBiomodProjFuture50,
    selected.models = models,
	binary.meth=metrics,
	filtered.meth=metrics,
	compress=TRUE)
  

  do.call(file.remove,list(list.files(pattern="temp*"))) 
}

save.image()


#myModelsOutlapply <-lapply(sp.n,BioModApply)

## alternatively with lapply
# myLapply_SFModelsOut <- lapply( sp.n, BioModApply)

# sph.umich.edu/biostat/computing/cluster/examples/r.html
library(Rmpi)

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
#library(snow)
# https://www.osc.edu/~kmanalo/r_parallel
#slaves <- as.numeric(Sys.getenv(c("PBS_NP")))-1
#slaves <- mpi.universe.size() - 1

#cl <- makeCluster(slaves, type = "MPI")
# or following http://homepage.divms.uiowa.edu/~luke/R/cluster/cluster.html
#cl <- getMPIcluster()
#clusterExport(cl, SLURMdat)

#tick <- proc.time()
#mySFModelsOut<-parLapplyLB(cl,sp.n,BioModApply)
#tock <- proc.time() - tick

#cat("\nsnow w/ Rmpi test times using", slaves, "MPI slaves: \n")

#stopCluster(cl)
#mpi.quit()

########

# snow-modelling.R
library(Rmpi)
library(snow)

# via https://rcc.uchicago.edu/docs/software/environments/R/index.html#multicore
np <- mpi.universe.size() - 1
cluster <- makeMPIcluster(np)
ClusterExport(cl=cluster, list("sp.n", "BioModApply","presmodstack","f50modstack","f70modstack","metrics","models"),
              envir=environment())

ClusterApply(cl, library(biomod2,snow,Rmpi))
applyOut<-parLapply(cluster, sp.n,BioModApply)
save(applyOut)
stopCluster(cluster)
mpi.exit()


########
#cl <- makeCluster(64, type = "MPI")
#clusterExport(cl, SLURMdat)

q
#mySFModelsOut<-clusterApply(cl,sp.n,fun=BioModApply)
#mySFModelsOut<-parLapplyLB(cl,sp.n,BioModApply)
#mySFModelsOut<-parLapply(cl,sp.n,BioModApply)
#stopCluster(cl)
#Rmpi:mpi.finalize()
#######




#######
# snowfall initialization
# library(snowfall)
# cl <- makeCluster(length(sp.n), type = "MPI") 
# sfInit(parallel=TRUE,cpus=27)

## Export packages to snowfall
# sfLibrary('biomod2', character.only=TRUE)
# sfExportAll()

# you may also use sfExportAll() to export all your workspace variables
## Do the run
# mySFModelsOut <- sfLapply(sp.n,BioModApply)

## stop snowfall

# sfStop( nostop=FALSE )
#######

