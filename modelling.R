
# establish directories
root<-"E:\\thesis"
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

# 
dir_clim<-paste0(dir_dat,"/clim")
dir_pres<-paste0(dir_clim,"/present")
dir_fut<-paste0(dir_clim,"/future")

dir_p.mosaics<-paste0(dir_pres,"/2.0/")
dir_f.mosaics<-paste0(dir_fut,"/1.4/")

dir_stacks<-paste0(dir_clim,"/stacks/")

## recursively create directories if not already there
folders<-as.list(ls())
for (i in folders[[i]])  { 
  folder<-get(i)
  dir.create(folder,recursive=TRUE,pattern="dir") 
} 


## read in functions
# Function to Install and Load R Packages
# borrowed from Pratik Patil, https://stackoverflow.com/questions/9341635/check-for-installed-packages-before-running-install-packages
Install_And_Load <- function(Required_Packages)
{
  Remaining_Packages <- Required_Packages[!(Required_Packages %in% installed.packages()[,"Package"])];
  
  if(length(Remaining_Packages)) 
  {
    install.packages(Remaining_Packages);
  }
  for(package_name in Required_Packages)
  {
    library(package_name,character.only=TRUE,quietly=TRUE);
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
                     "snowfall"))

Install_And_Load(requiredPackages)

#####################################
#       BIOMOD formatting Data      #
#####################################

# set wd inside biomod subdirectory
setwd(dir_bm)


# comment out when running cluster
pa<-read.csv(file=paste0(dir_out,"/pa_dataframe.csv"))
pa<-data.frame(pa)


# get vector of species names
names<-paste0(colnames(pa))
sp.n= dput(names[49]) #vector of species name(s), excluding lat and long cols

#sp.n=c("TuxpeÃ±o"
#,"Arrocillo.Amarillo"
#        )     

presmodstack<-stack(paste0(dir_stacks,"present_modstack.grd"),paste0(dir_dat,"/ind/pob-ind.grd"))
f50modstack<-stack(paste0(dir_stacks,"f50_modstack.grd"),paste0(dir_dat,"/ind/pob-ind.grd"))
f70modstack<-stack(paste0(dir_stacks,"f70_modstack.grd"),paste0(dir_dat,"/ind/pob-ind.grd"))
# load(paste0(dir_bm,"/.RData"))

library(quickPlot)
names(presmodstack)<-gsub("crop_wc2.0_","",layerNames(presmodstack))
names(presmodstack)<-gsub("_30s","",layerNames(presmodstack))

names(f50modstack)<-layerNames(presmodstack)
names(f70modstack)<-layerNames(presmodstack)

setwd(dir_bm)


BioModApply <-function(sp.n) {
  rasterOptions()$tmpdir
  rasterOptions(tmpdir=root)
  maxentjar<-paste0(dir_R,"/maxent/maxent.jar")
  setwd(dir_bm)
  library(biomod2)
  
  library(mgcv)
  options(max.print=1000000)
  library(gbm)
  library(dismo)
  myRespName = sp.n
  myResp <- as.numeric(pa[,myRespName])
  myRespXY = pa[,c('Longitude.x.','Latitude.y.')]
  
  myExpl<-presmodstack
  myExplFuture50<-f50modstack
  myExplFuture70<-f70modstack
  
  
  # Barbet-Massin et al 2012:
  #   Overall, we recommend the use of a large number (e.g. 10 000) of pseudo-absences with equal
  # weighting for presences and absences when using regression techniques (e.g. generalised linear
  #                                                                        model and generalised additive model); averaging several runs (e.g. 10) with fewer pseudo-absences
  # (e.g. 100) with equal weighting for presences and absences with multiple adaptive regression splines
  # and discriminant analyses; and using the same number of pseudo-absences as available presences
  # (averaging several runs if few pseudo-absences) for classification techniques such as boosted regression
  # trees, classification trees and random forest. In addition, we recommend the random selection
  # of pseudo-absences when using regression techniques and the random selection of geographically
  # and environmentally stratified pseudo-absences when using classification and machine-learning
  # techniques
  
  
  # format input data for biomod
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       resp.xy = myRespXY,
                                       expl.var = myExpl,
                                       resp.name = myRespName,
                                       PA.nb.rep = 1,
                                       PA.nb.absences = 500,
                                       PA.strategy = "random",
                                       na.rm=TRUE
  )
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  # define model options
  
  
  # print default biomod options 
  default_ModelOptions <-BIOMOD_ModelingOptions()
  print(default_ModelOptions)
  
  install.library("gam")
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
  
  install.packages("ENMeval")
  library(caret)
  library(ENMeval)
  
  # download new version of code from Frank Breiner (writer of BIOMOD_tuning), attached here: http://r-forge.wu.ac.at/forum/forum.php?max_rows=75&style=nested&offset=152&forum_id=995&group_id=302
  
  source(paste0(dir_R,"/maices-enm/BIOMOD.tuning_v6.R"))
  library(doParallel);cl<-makeCluster(8);registerDoParallel(cl) 
  # devtools::install_github('topepo/caret/pkg/caret')
  library(caret)
  BIOMOD_TunedOptions <- BIOMOD_tuning(myBiomodData,
                                 env.ME = myExpl,
                                 n.bg.ME = ncell(myExpl)
                                 )
  stopCluster(cl)
  BIOMOD_ModelOptions<-Biomod.tuning$models.options
  

  capture.output(BIOMOD_TunedOptions$models.options,file=paste0(dir_out,"/model-opts/",myRespName,"_tuned_opts.txt"))
  
  # modeling
  myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData, 
    models = c("GLM","GAM","GBM","ANN","CTA","RF","MARS","FDA","MAXENT.Phillips",'MAXENT.Tsuruoka'), 
    models.options = Biomod.tuning$models.options, 
    NbRunEval=2,
    DataSplit=70,
    VarImport=2,
    models.eval.meth = c( 'KAPPA', 'TSS'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = TRUE,
    modeling.id = paste0(myRespName,"_current"))
  
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  print(paste0("Done Running Models for ",sp.n))
  
  dir.create(paste0(dir_out,"/",myRespName))
  
  
  # write data used for modelling
  capture.output(get_formal_data(myBiomodModelOut),
                 file=paste0(dir_out,"/",myRespName,"/",myRespName,"_model_data.txt"))
  ### eval current model
  
  print(paste0("Capturing Model Evaluations for ",myRespName))
  evalmods<-get_evaluations(myBiomodModelOut,as.data.frame=TRUE)
  write.csv(evalmods,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_models_eval.csv"))
  
  ### get variable importance
  modevalimport<-get_variables_importance(myBiomodModelOut,as.data.frame=TRUE)
  write.csv(modevalimport,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_var_imp.csv"))
  
  ### get model summaries
  capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_current/",myRespName,"_PA1_Full_ANN",sep="")))))
                 ,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_ANN_summary.txt"))
  
  capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_current/",myRespName,"_PA1_Full_CTA",sep="")))))
                 ,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_CTA_summary.txt"))
  
  # fda1<-get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"_current/",myRespName,"_PA1_Full_FDA",sep=""))))
  # capture.output(as.table(fda1$confusion),file==paste0(dir_out,"/",myRespName,"/",myRespName,"_FDA-conf_summary.csv"))

    capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"_current/",myRespName,"_PA1_Full_GAM",sep="")))))
                 ,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_GAM_summary.txt"))
    
  capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"_current/",myRespName,"_PA1_Full_GBM",sep="")))))
                 ,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_GBM_summary.txt"))
  
  capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"_current/",myRespName,"_PA1_Full_GLM",sep="")))))
                 ,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_GLM_summary.txt"))
  
  capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"_current/",myRespName,"_PA1_Full_MARS",sep="")))))
                 ,file=paste0(dir_out,"/",myRespName,"/",myRespName,"_MARS_summary.txt"))
  
  copy()
  # summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_PA1_Full_MAXENT.Phillips",sep="")))))
  # maxent_t1<-get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_PA1_Full_MAXENT.Tsuruoka",sep=""))))
  
  models<-c("GLM","GAM","GBM","ANN","CTA","RF","MARS","FDA",'MAXENT.Tsuruoka')
  for (mods in models){
    loadedmodel<-biomod2::BIOMOD_LoadModels(myBiomodModelOut,models=mods)
    dir.create(paste0(dir_figs,"/",myRespName))
    response.plot2(models= loadedmodel,
                   Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                   show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
                   save.file = "tiff",
                   name=paste0(dir_figs,"/",myRespName,"/",myRespName,"_",mods,"_curves"),
                   col = c("blue", "red"),
                   legend = TRUE,
                   data_species = get_formal_data(myBiomodModelOut,'resp.var'),
                   ImageSize=1000
    )
    
  }

  
  # get model scores by eval metrics
  dir.create(paste0(dir_figs,"/",myRespName))
  png(filename=paste0(dir_figs,"/",myRespName,"/",myRespName,"_model_scores-kappa-tss.png"))
  models_scores_graph(myBiomodModelOut,
                      metrics = c( 'KAPPA', 'TSS'),
                      by = 'models',
                      plot = TRUE)
  dev.off()
  
  png(filename=paste0(dir_figs,"/",myRespName,"/",myRespName,"_observations.png"))
  presab<-myResp
  presab[is.na(presab)]<-0
  level.plot(data.in = presab,
             XY = myRespXY )
  dev.off()
  
  # model projections
  myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExpl,
    proj.name = 'current',
    selected.models = 'all',
    binary.meth = c( 'KAPPA', 'TSS'),
    compress = TRUE,
    clamping.mask = TRUE,
    output.format = '.grd')
  
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  
  
  coutputFolderName <- "proj_current"
  coutputFolder <- paste(myRespName,coutputFolderName,sep="/")      #FINDING DIRECTORY FOR ENSEMBLE MODEL OUTPUT FOR REPORT GENERATION
  print("Done Projecting Models")
  
  
  # mod_proj <- get_predictions(myBiomodProj) 
  # plot(mod_proj)
  # plot(subset(mod_proj,2), main = "RF projections") # select layer in stack to plot
  
  # ensemble modeling
  myBiomodEM <- BIOMOD_EnsembleModeling(
    modeling.output = myBiomodModelOut,
    chosen.models = 'all',
    em.by="PA_dataset+repet",
    eval.metric = c( 'KAPPA', 'TSS'),
    eval.metric.quality.threshold = c(rep(0.6,2)),
    prob.mean = T,
    prob.cv = T,
    prob.ci = T,
    prob.ci.alpha = 0.05,
    prob.median = T,
    committee.averaging = T,
    prob.mean.weight = T,
    prob.mean.weight.decay = 'proportional' )
  
  
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  ### eval current model
  
  print(paste0("Capturing Model Ensemble Evaluations for ",sp.n))
  enevalmods<-get_evaluations(myBiomodEM,as.data.frame=TRUE)
  write.csv(enevalmods,file=paste0(dir_out,"/eval/formal_ensemble_evaluation.csv"))
  
  
  # current ensemble projection
  myBiomodEF <- BIOMOD_EnsembleForecasting(
    EM.output = myBiomodEM,
    projection.output = myBiomodProj,
    binary.meth = c( 'KAPPA', 'TSS'),
    )
  
  myBiomodEF
  # EF_stack<- raster::stack(paste0(dir_bm,"/",sp.n,"/proj_current/proj_current_",sp.n,"_ensemble.grd"))
  # plot a layer in the ensemble stack
  # plot(EF_stack[[6]])
  
  
  # future projections for rcp 85 period 70
  myBiomodProjFuture70 <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExplFuture70,
    proj.name = 'rcp85_70',
    selected.models = 'all',
    binary.meth = c( 'KAPPA', 'TSS'),
    compress = 'xz',
    clamping.mask = T,
    output.format = '.grd')
  
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  
  # future projections for rcp 85 period 50
  myBiomodProjFuture50 <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExplFuture50,
    proj.name = 'rcp85_50',
    selected.models = 'all',
    binary.meth = c( 'KAPPA', 'TSS'),
    compress = 'xz',
    clamping.mask = T,
    output.format = '.grd')
  
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  f70BiomodEF <- BIOMOD_EnsembleForecasting(
    EM.output = myBiomodEM,
    projection.output = myBiomodProjFuture70)
  cat("\n\nExporting Ensemble as grd ...\n\n")
  
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  f50BiomodEF <- BIOMOD_EnsembleForecasting(
    EM.output = myBiomodEM,
    projection.output = myBiomodProjFuture50)
  cat("\n\nExporting Ensemble as grd ...\n\n")
  
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  # myBiomodEF
  # EF_70stack<- raster::stack(paste0(dir_bm,"/",sp.n,"/proj_current/proj_rcp85_70",sp.n,"_ensemble.grd"))
  # plot a layer in the ensemble stack
  # plot(EF_stack[[6]])
  
  # myBiomodEF
  # EF_50stack<- raster::stack(paste0(dir_bm,"/",sp.n,"/proj_current/proj_rcp85_50",sp.n,"_ensemble.grd"))
  # plot a layer in the ensemble stack
  # plot(EF_stack[[6]])
  
  projs<-c("proj_current","proj_rcp85_50","proj_rcp85_70")
  
  for (p in projs){
    
  
  ensemblestack<-raster::stack(paste0(dir_bm,"/",sp.n,"/",p,"/",p,"_",sp.n,"_ensemble.grd"))
  ensembleconsensus<-stackApply(ensemblestack,indices=rep(1:1,nlayers(ensemblestack)),fun=mean)
  ensembleconsensus100<-ensembleconsensus/1000
  writeRaster(ensembleconsensus100,file=paste0(dir_bm,"/",sp.n,"/",p,"/",p,"_",sp.n,"_mean_consensus.grd"),format="raster",overwrite=T)
  }
 
  
#EXPORTING ENSEMBLE MODEL PROJECTION AS ASCII FOR USE IN OUTSIDE MAPPING SOFTWARE
gridName = paste(coutputFolderName,myRespName,"ensemble.grd",sep="_")  
gridDir = paste(myRespName,coutputFolderName,gridName,sep="/")
MyRaster = raster(paste(myRespName,coutputFolderName,gridName,sep="/"))
writeRaster(MyRaster,file=paste0(sp.n,"_EnsembleRaster.asc"), format = 'ascii', overwrite = TRUE)
print("Done Exporting Ensemble as ASCII")

cat("\n\nExporting Model Plots ...\n\n")
#EXPORTING PLOTS FOR EACH MODEL & ENSEMBLE TO 'PLOTS' FOLDER IN WD

model_names <- c("GLM","GBM","GAM","ANN","SRE","CTA","RF","MARS","FDA","MAXENT.Phillips",'MAXENT.Tsuruoka')

              
model_names <- c("GLM","GBM","GAM","ANN","SRE","CTA","RF","MARS","FDA","MAXENT.Phillips",'MAXENT.Tsuruoka')
i<-model_names[1]
for (i in myBiomodProj@models.projected){
  raster::plot(i)
}

png(paste0("plots/ensemble_",sp.n,".png"))
raster::plot(myBiomodEF)
dev.off()
print("Done Exporting Model Plots")

# }

# define a mask of studied
alphaMap <- reclassify(subset(myExpl,1), c(-Inf,Inf,0))
currentensemble<-stack(file.path(variety,
                                 "proj_current", 
                                 paste("proj_current_",
                                       variety, 
                                       "_ensemble.grd", sep="")))
names(currentensemble)

# # add all other species map
for(variety in sp.n){
  # add layer
  alphaMap <- 
    alphaMap + 
    subset(stack(file.path(variety,
                           "proj_current", 
                           paste("proj_current_",
                                 variety, 
                                 "_TSSbin.grd", sep=""))), 1)
}

# summary of created raster
plot(alphaMap)


