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
library(mgcv)


#################################################################
# PRELIMINARY DATA FORMATTING
#################################################################

print("Reading in Presence/Absence Matrix")

### read in PAM with projected coordinates
pa<-data.frame(read.csv(file=paste0(dir_out,"/pa_dataframe_EA.csv")))
colnames(pa) 

# filter species by obs count
print("Filtering varieties with > 14 observations:")
i <- colSums(pa[,4:ncol(pa)],na.rm=T) > 14 
pa<-pa[,i]
colnames(pa)

# drop unwanted columns by name to get vector of species names
#colnames(pa)
drops <- c("X","Longitude.x.","Latitude.y.","X_Lambert", "Y_Lambert")
pa_sp.n<-pa[ , !(names(pa) %in% drops)]

colnames(pa_sp.n)


names<-paste0(colnames(pa_sp.n))
sp.n= dput(names) # keep only species name, remove lat/long/etc. 

print("The species to be modeled:")
print("")
sp.n

###################################
## original version, untransformed
#i <- (colSums(pa[4:ncol(pa)],na.rm=T)) > 14 # filter species by obs count
#pa<-pa[,i]
#nspec<-ncol(pa[,4:ncol(pa)]) # number of species modelled
#print(paste0(nspec," species/varieties selected in PA matrix"))

#names<-paste0(colnames(pa))
#names

#sp.n= dput(names [c(4:length(names))] # keep only species name, remove lat/long/etc. 
#) #vector of species name(s), excluding lat and long cols
#Encoding(sp.n)<-"latin1"
#sp.n= sp.n[1]

#sp.n= dput(names [c(4:11)] # keep only species name, remove lat/long/etc. 
#) #vector of species name(s), excluding lat and long cols
###################################

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

presmodstack<-stack(paste0(dir_stacks,"present_modstack_EA1km.grd"))

f50modstack<-stack(paste0(dir_stacks,"f50_modstack_EA1km.grd"))
f70modstack<-stack(paste0(dir_stacks,"f70_modstack_EA1km.grd"))

names(presmodstack)

print("Reading in estimated prevalence .csv:")

prev<-as.data.frame(read.csv(file.path(dir_maices,"tau.csv")))
colnames(prev)
prev$Raza_prima
prev$Raza_prima <- sub(" ", "_", prev$Raza_prima)
prev$Raza_prima <- sub(" ", "_", prev$Raza_prima)
prev$Raza_prima <- sub(" ", "_", prev$Raza_prima)
prev$Raza_prima

#################################################################
# INITIALIZE FUNCTION TO APPLY TO EACH VARIETY
#################################################################


# following https://rpubs.com/dgeorges/190889
print("Setting Working Directory:")
setwd(dir_bm) 

#maxent.background.dat.dir <- paste0(getwd(),"/maxent_bg")
#dir.create(maxent.background.dat.dir, showWarnings = FALSE, recursive = TRUE)

## resave explanatory data

#for(var_ in names(presmodstack)){

#	cat("\n> saving", paste0(var_, ".asc"))
#
#	f<-file.path(maxent.background.dat.dir, paste0(var_, ".asc"))
#	#if (!file.exists(f))

#		writeRaster(subset(presmodstack, var_), 
#               	filename = f,
#              	overwrite = TRUE)
#}



#options(max.print=1000000)  # set max.print option high to capture outputs
#path.to.maxent.jar<-file.path(getwd(),"maxent.jar") # define maxent jar location


print("Establishing Function BioModApply() :")


#################################################################
# DEFINE MODELING FUNCTION
#################################################################
BioModApply <-function(sp.n) {
  tryCatch({
    setwd(dir_bm)

    # keep temporary files separate for sp.n 
    #... the maxent output and raster tmp files might get overwritten and cause errors in parallel


    #rasterOptions()$tmpdir      # get raster temp director
    rtmpdir<-paste0(root,"/tmp/",sp.n)
    dir.create(rtmpdir, showWarnings = FALSE, recursive = TRUE)
    rasterOptions(tmpdir=rtmpdir)  # set raster temp directory
    tempfile(pattern=sp.n) ## pattern for temp files prefixed with sp.n
	
    ### 
    # if using raster stack
    myRespName = sp.n
    #myRespName = sp.n[1]
    myResp <- as.numeric(pa[,myRespName])
    myRespXY = pa[,c('Longitude.x.','Latitude.y.')] # actually x, y in lambert conformal # letsR::presab points makes column names long/lat
    myExpl<-presmodstack
    myExpl_proj<-presmodstack
    myExplFuture50<-f50modstack
    myExplFuture70<-f70modstack
    
    
    
    ###
    
    ###
    # if using matrix
    #myRespName = sp.n
    # myRespName = sp.n[1]
    #myResp <- as.numeric(pa[,myRespName])
    #myRespXY = pa[,c('Longitude.x.','Latitude.y.')] # actually x, y in lambert conformal # letsR::presab points makes column names long/lat
    #myExpl<-extract(presmodstack,as.data.frame(myRespXY),na.rm=FALSE)
    #myExpl_proj<-raster::as.matrix(presmodstack)
    #myExplFuture50<-raster::as.matrix(f50modstack)
    #myExplFuture70<-raster::as.matrix(f70modstack)
    

    ###

    # if using matrix as predictor
    allmodels<-c("GLM","GAM","GBM","ANN","CTA","RF","MARS","FDA","MAXENT.Phillips","MAXENT.Tsuruoka")
    models = c("MAXENT.Phillips")
    metrics = c(  'KAPPA', 'TSS', 'ROC', 'FAR','SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS')
    
    #load("/gpfs/home/scg67/thesis/02_R/maices-enm/.RData")

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
                                         resp.xy = myRespXY, # if using matrix comment this line out
                                         expl.var = myExpl,
                                         resp.name = myRespName,
                                         na.rm=TRUE,
                                         PA.nb.rep =1,
                                         PA.nb.absences = 20000,
                                         PA.strategy = 'random'
    )

    #do.call(file.remove,list(list.files(pattern="temp*"))) 
    
    #################################################################
    # DEFINE MODEL OPTIONS
    #################################################################
    
    
    # print default biomod options 
    #default_ModelOptions <-BIOMOD_ModelingOptions()
    #print(default_ModelOptions)
    
    library(gam)

    # edit default options accordingly
    BIOMOD_ModelOptions <-
      BIOMOD_ModelingOptions(
        
        GLM = list(
          type = 'quadratic',
          interaction.level = 0,
          myFormula = NULL,
          test = 'AIC',
          # tk BIC?
          family = binomial(link = 'logit'),
          mustart = 0.5,
          control = glm.control(
            epsilon = 1e-08,
            maxit = 50,
            trace = FALSE
          )
        ),
        
        
        GBM = list(
          distribution = 'bernoulli',
          n.trees = 2500,
          interaction.depth = 7,
          n.minobsinnode = 5,
          shrinkage = 0.001,
          bag.fraction = 0.5,
          train.fraction = 1,
          cv.folds = 3,
          keep.data = FALSE,
          verbose = FALSE,
          perf.method = 'cv'
        ),
        
        CTA = list(
          method = 'class',
          parms = 'default',
          cost = NULL,
          control = list(
            xval = 5,
            minbucket = 5,
            minsplit = 5,
            cp = 0.01,
            maxdepth = 25
          )
        ),
        
        
        ANN = list(
          NbCV = 5,
          size = NULL,
          decay = NULL,
          rang = 0.1,
          maxit = 200
        ),
        
        SRE = list(quant = 0.025),
        
        FDA = list(method = 'mars',
                   add_args = NULL),
        
        MARS = list(
          type = 'simple',
          interaction.level = 0,
          myFormula = NULL,
          nk = NULL,
          penalty = 2,
          thresh = 0.001,
          nprune = NULL,
          pmethod = 'backward'
        ),
        
        RF = list(
          do.classif = TRUE,
          ntree = 500,
          mtry = 'default',
          nodesize = 5,
          maxnodes = NULL
        ),
        
        MAXENT.Phillips = list(
          #memory_allocated = 10240,
          #background_data_dir = maxent.background.dat.dir, # https://rpubs.com/dgeorges/190889
          #maximumbackground = 10000,
          maximumiterations = 10000,
          visible = FALSE,
          linear = FALSE,
          quadratic = FALSE,
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
          defaultprevalence = 0.1
        ),
        
        MAXENT.Tsuruoka = list(
          l1_regularizer = 0,
          l2_regularizer = 0,
          use_sgd = FALSE,
          set_heldout = 0,
          verbose = FALSE
        )
      )
    
    #bm.opt.maxent.bg <- BIOMOD_ModelingOptions(MAXENT.Phillips = list(#memory_allocated = 8192,
    #background_data_dir = maxent.background.dat.dir, # https://rpubs.com/dgeorges/190889
    #maximumbackground = 10000,
    #maximumiterations = 5000,
    #visible = FALSE,
    #linear = FALSE,
    #quadratic = TRUE,
    #product = FALSE,
    #threshold = FALSE,
    #hinge = TRUE,
    #lq2lqptthreshold = 80,
    #l2lqthreshold = 10,
    #hingethreshold = 10,
    #beta_threshold = -1,
    #beta_categorical = -1,
    #beta_lqp = -1,
    #beta_hinge = 0.5,
    #betamultiplier = 1,
    #defaultprevalence = 0.5,
    #path_to_maxent.jar = path.to.maxent.jar))
    
    #################################################################
    # TUNE MODEL OPTIONS -- NOT WORKING RN
    #################################################################
    
    # install.packages("ENMeval")
    #library(doParallel);cl<-makeCluster(8);registerDoParallel(cl) 
    # devtools::install_github('topepo/caret/pkg/caret')
    #library(caret)
    # download new version of code from Frank Breiner (author), attached here: http://r-forge.wu.ac.at/forum/forum.php?max_rows=75&style=nested&offset=152&forum_id=995&group_id=302
    
    #library(dismo) #ensure that maxent.jar is also in R/lib/dismo/java
    #library(ENMeval)
    
    #source(paste0(dir_R,"/maices-enm/BIOMOD.tuning_v6.R"))
    #BIOMOD_TunedOptions <- BIOMOD_tuning(
    #  myBiomodData,
    #  models = "MAXENT.Phillips",
    #  env.ME = myExpl,
    #  models.options = bm.opt.maxent.bg,
    #  n.bg.ME = ncell(myExpl),
    #  metric.ME = "ROC"
    #)
    
    #if (exists("BIOMOD_TunedOptions"))
    #{
    #  BIOMOD_TunedOptions <- get("BIOMOD_TunedOptions")
    #  BIOMOD_ModelOptions <- BIOMOD_TunedOptions$models.options
    #  
    #}
    
    
    #capture.output(BIOMOD_TunedOptions$models.options,file=paste0(dir_out,"/model-opts/",myRespName,"_tuned_opts.txt"))
    
    #################################################################
    # BUILD MODELS
    #################################################################
    
    ### 
    # https://r-forge.r-project.org/forum/message.php?msg_id=41945&group_id=302
    # 
    # "The way you are using the Yweight arg is the right one, but not all biomod2 models suppor currently weights arguments and MAXENT only use weights via prevalence argument (not directly weights given via Yweights).
    # Here a summary of models that use or not Yweight arg :
    
    #===========================
    # model weigths supported ?
    # GLM => YES
    # GBM => YES
    # GAM => YES
    # CTA => YES
    # ANN => YES
    # SRE => NO
    # FDA => YES
    # MARS => YES
    # RF => NO
    # MAXENT => NO (weights only via prevalence arg) "
    
    ####
    # this would have been nice though,
    # following: https://www.researchgate.net/post/Can_I_run_a_raster_bias_layer_in_biomod22
    
    ##import raster(bias_layer)]
    #myBiasfile
    #myBiasfile<-raster(paste0(dir_stacks,"/bias.grd"))
    
    #myPresPAdf <- data.frame(myBiomodData@coord,obs = myBiomodData@data.species, myBiomodData@PA)
    #head(myPresPAdf,100)
    #summary(myBiomodData)
    ## add a random weight vector
    #myPresPAdf$yweights <- extract(bias,myBiomodData@coord)
    #range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    
    #myPresPAdf$yweights_stnd<-range01(myPresPAdf$yweights)
    #head(myPresPAdf)
    #sp_weights<-myPresPAdf$yweights_stnd
    
    
    
    myBiomodModelOut <-
      BIOMOD_Modeling(
        myBiomodData,
        models = c("MAXENT.Phillips"),
        models.options = BIOMOD_ModelOptions,
        NbRunEval = 25,
        DataSplit = 70,
        Prevalence = prev$tau[prev$Raza_prima == myRespName],
        # use estimated tau using Perales et al 2015 outputs 
        VarImport = 5,
        models.eval.meth = metrics,
        SaveObj = TRUE,
        rescal.all.models = FALSE,
        do.full.models = FALSE,
        modeling.id = paste0(myRespName, "_current")
      )
    
    #do.call(file.remove,list(list.files(pattern="temp*"))) 
    
    #print(paste0("Done Running Models for ",sp.n))
    

    if (exists("myBiomodModelOut"))
    {
      myBiomodModelOut <- get("myBiomodModelOut")
      
    
    #################################################################
    # CAPTURE DATA INPUT AND MODEL OUTPUTS
    #################################################################
    
	dir.create(paste0(dir_out, "/", myRespName),showWarnings = FALSE)
    
    # write data used for modelling
    capture.output(
      get_formal_data(myBiomodModelOut),
      file = paste0(getwd(),"/",myRespName,"/",myRespName, "_model_data.txt")
    )
    
    # get model evaluations by metrics
    evalmods <- get_evaluations(myBiomodModelOut, as.data.frame = TRUE)
    evalmods$variety <- myRespName
    write.csv(evalmods,
              file = paste0(getwd(),"/",myRespName,"/",myRespName, "_models_eval.csv"))
    
    ### get variable importance
    modevalimport <-
      get_variables_importance(myBiomodModelOut, as.data.frame = TRUE)
    colnames(modevalimport)
    
    write.csv(modevalimport,
              file = paste0(getwd(),"/",myRespName,"/",myRespName, "_var_imp.csv"))
    
    
	
    # #################################################################
    # BUILD ENSEMBLE MODELS
    #################################################################
    #print(paste0("Building Ensemble Models for ",myRespName))
    
    # ensemble modeling
    myBiomodEM <- BIOMOD_EnsembleModeling(
      modeling.output = myBiomodModelOut,
      chosen.models = 'all',
      em.by = 'all',
      eval.metric = metrics,
      eval.metric.quality.threshold = NULL,
      prob.mean = F,
      prob.cv = T,
      prob.ci = T,
      prob.ci.alpha = 0.05,
      prob.median = T,
      committee.averaging = T,
      prob.mean.weight = T,
      prob.mean.weight.decay = "proportional",
      VarImport = 10
    )
    
    #################################################################
    # CAPTURE ENSEMBLE MODEL OUTPUTS
    #################################################################
    
    #print(paste0("Capturing Ensemble Model Outputs for  ",myRespName))
    
    # write em models built
    capture.output(
      get_built_models (myBiomodEM),
      file = paste0(myRespName, "_em_models.txt")
    )
    
    # capture em model evals
    #print(paste0("Capturing Ensemble Models Evaluations ",myRespName))
    evalmodsEM <- get_evaluations(myBiomodEM, as.data.frame = TRUE)
    evalmodsEM$variety <- myRespName
    write.csv(evalmodsEM,
              file = paste0(getwd(),"/",myRespName,"/",myRespName,"_models_EM_eval.csv"))
    
    #do.call(file.remove,list(list.files(pattern="temp*")))
    
    ### eval current model
    
    #print(paste0("Capturing Model Ensemble Evaluations for ",myRespName))
    enevalmods <- get_evaluations(myBiomodEM, as.data.frame = TRUE)
    write.csv(enevalmods,
              file = paste0(getwd(),"/",myRespName,"/",myRespName, "_em_evals-df.csv"))
    
    } # if.exists biomodout > ensemble, project, then forecast
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n", sep=" ", master=TRUE )})
  
}
 

#################################################################
# PROJECTIONS AND FORECASTING
#################################################################
BioProjectApply <-function(sp.n) {
  tryCatch({
    setwd(dir_bm)
    
    setwd(dir_bm) 
    
    myExpl<-presmodstack
    myExplFuture50<-f50modstack
    myExplFuture70<-f70modstack
    
    
    bm_out_file <- load(paste0(getwd(),"/",sp.n,"/",sp.n,".",sp.n,"_current.models.out"))
    myBiomodModelOut <- get(bm_out_file)

    # model projections
    myBiomodProj <- BIOMOD_Projection(
      modeling.output = myBiomodModelOut,
      new.env = myExpl_proj,
      proj.name = 'current',
      selected.models = 'all',
      binary.meth = metrics,
      compress = "gzip",
      build.clamping.mask = T,
      output.format = '.RData',
      do.stack = FALSE,
      keep.in.memory = FALSE,
      silent = TRUE
    )
    
    #print(paste0("Projecting  onto Future (2050) for ",myRespName))
    # future projections for rcp 85 period 50
    myBiomodProjFuture50 <- BIOMOD_Projection(
      modeling.output = myBiomodModelOut,
      new.env = myExplFuture50,
      proj.name = 'rcp85_50',
      selected.models = 'all',
      binary.meth = metrics,
      compress = "gzip",
      build.clamping.mask = T,
      output.format = '.RData',
      do.stack = FALSE,
      keep.in.memory = FALSE,
      silent = TRUE
    )
    
    #print(paste0("Projecting onto Future (2070) for ",myRespName))
    # future projections for rcp 85 period 70
    myBiomodProjFuture70 <- BIOMOD_Projection(
      modeling.output = myBiomodModelOut,
      new.env = myExplFuture70,
      proj.name = 'rcp85_70',
      selected.models = 'all',
      binary.meth = metrics,
      compress = "gzip",
      build.clamping.mask = T,
      output.format = '.RData',
      do.stack = FALSE,
      keep.in.memory = FALSE,
      silent = TRUE
    )
    
    #unlink(rtmpdir,recursive=TRUE)
    
    #print(paste0("Performing Ensemble Forcasting onto Current Data for ",myRespName))
    #################################################################
    # GET ENSEMBLE FORECASTS
    #################################################################
    

    bm_emout_file <- load(paste0(getwd(),"/",sp.n,"/",sp.n,".",sp.n,"_currentensemble.models.out"))
    myBiomodEM <- get(bm_emout_file)

    
    # current ensemble projection
    myBiomodEF <- BIOMOD_EnsembleForecasting(
      EM.output = myBiomodEM,
      projection.output = myBiomodProj,
      selected.models = 'all',
      binary.meth = metrics,
      compress = "gzip",
      total.consensus = TRUE
    )
    
    #print(paste0("Performing Ensemble Forcasting onto Future (2070) Data for ",myRespName))
    
    f70BiomodEF <- BIOMOD_EnsembleForecasting(
      EM.output = myBiomodEM,
      projection.output = myBiomodProjFuture70,
      selected.models = 'all',
      binary.meth = metrics,
      compress = "gzip",
      total.consensus = TRUE
    )
    
    #cat("\n\nExporting Ensemble as grd ...\n\n")
    
    #do.call(file.remove,list(list.files(pattern="temp*")))
    
    #print(paste0("Performing Ensemble Forcasting onto Future (2050) Data for ",myRespName))
    
    f50BiomodEF <- BIOMOD_EnsembleForecasting(
      EM.output = myBiomodEM,
      projection.output = myBiomodProjFuture50,
      selected.models = 'all',
      binary.meth = metrics,
      compress = "gzip",
      total.consensus = TRUE
    )
            
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n", sep=" ", master=TRUE )})
  
}


#######
# snowfall initialization
#install.packages("Rmpi", lib="/gpfs/home/scg67/R/x86_64-pc-linux-gnu-library/3.4/", repos='http://cran.us.r-project.org')
library(snowfall)
library(Rmpi)
library(raster)
#args = commandArgs(trailingOnly = TRUE);
#ncpus = args[1];

# ------------------------------------------------------------------------
# initialize parallel mode using SLURM and command-line args
# ------------------------------------------------------------------------

# # get number of SLURM tasks
ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS"))

print(ntasks)
print(class(ntasks))
print(paste0(ntasks, " tasks"))

# # get number of CPUs per node from SLURM
ncputask <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
print(paste0(ncputask, " N CPUs per task"))

# # total cpus equals tasks * cpus per taks
ncpus <- ncputask * ntasks
print(paste0(ncpus, " CPUs"))

# other options:
#ncpus<-detectCores()
#ncpus<-future::availableCores()
#ncpus<-mpi.universe.size()-1

sfSetMaxCPUs(ncpus)
# # initialize cluster
sfInit(parallel = TRUE,
       cpus = ncpus,
       type = "MPI",
       slaveOutfile = paste0(Sys.getenv("SLURM_JOB_ID"),"_sfInit.log"))

# # Export packages to snowfall
print("Exporting packages to cluster")
sfLibrary('biomod2', character.only = TRUE)
sfLibrary('raster', character.only = TRUE)


# # Export environment to snowfall
print("Exporting environment to cluster")
#sfExportAll()

sfExport("sp.n","BioModApply","presmodstack","f50modstack","f70modstack",
         "pa","prev","dir_out","dir_bm","root")

# # Do the run
# use sfClusterApplyLB when machines have different specs
print("Running >> sfLapply(sp.n,BioModApply)") # sfLapply
mySFModelsOut <- sfLapply(sp.n[1:length(sp.n)], BioModApply)

#mySFModelsOut <- sfLapply(sp.n[1:20], BioModApply)
#mySFModelsOut <- sfLapply(sp.n[21:40], BioModApply)
#mySFModelsOut <- sfLapply(sp.n[41:length(sp.n)], BioModApply)

# # Stop the cluster
print("Stopping cluster")

sfStop()
