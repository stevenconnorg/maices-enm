
# establish directories
root<-"~/thesis"
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

dir_stacks<-paste0(dir_clim,"/stacks/")

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
                     "snowfall",
                     "quickPlot",
                     "biomod2",
                     "mgcv",
                     "gbm",
                     "dismo"
))

library(biomod2)
library(snowfall)
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

eval_metrics<-c( 'KAPPA', 'TSS')
allmodels<-c("GLM","GAM","GBM","ANN","CTA","RF","MARS","FDA",'MAXENT.Phillips','MAXENT.Tsuruoka')



projects<-c("proj_current","proj_rcp85_50","proj_rcp85_70")




#################################################################
# LOADING MODELS AND PROJECTIONS AFTER MODELLING
#################################################################
### 
projects<-c("proj_current","proj_rcp85_50","proj_rcp85_70")


### UNDER CONSTRUCTION...

# for (p in projects){...}

# for (var in sp.n) {
# modelfile<-load(paste0(dir_bm,"/",var,"/",var,".",var,"_current.models.out"))
# ...


###...

# }
# 


# Pep.currentmods<-load(paste0(dir_bm,"/Pepitilla/Pepitilla.Pepitilla_current.models.out"))
# Pepitilla.current.mods<-get(Pep.currentmods)

# data <- cbind(Pepitilla=get_formal_data(Pepitilla.current.mods,'resp.var'), 
#               get_formal_data(Pepitilla.current.mods,'expl.var'))
# View(data)

# models<-BIOMOD_LoadModels(Pepitilla.current.mods)
# variables_importance(models,data=data)
sp.n[]


### READ IN PLOTS FROM FILE
# plot models
# current_mods<-load(paste0(dir_bm,"/",sp.n,"/",sp.n,".",sp.n,"_current.models.out"))
# current_mod<-get(current_mods)

# proj<-load(paste0(dir_bm,"/",sp.n,"/",sp.n,".",sp.n,"_current.models.out"))
# current_mod<-get(current_mods)


# get observation data formatted
pa<-read.csv(file=paste0(dir_out,"/pa_dataframe.csv"))
pa<-data.frame(pa)

names<-paste0(colnames(pa))

sp.n= dput(names   [66]) #vector of species name(s), excluding lat and long cols
projects<-c("proj_current","proj_rcp85_50","proj_rcp85_70")
models = c("GLM","GAM","GBM","ANN","CTA","RF","MARS","FDA","MAXENT.Phillips",'MAXENT.Tsuruoka')


# get individual projections
for (sp in sp.n){
  for (p in projects){
    
    sp<-sp.n
    p<-projects[1]
    # get individual projections
    for (mod in models){
      individualproj<-list.files(paste0(dir_bm,"/",sp,"/",p,"/individual_projections/"),pattern=".grd",full.names = T)
      for (proj in individualproj){
        raster(proj)
      }
      emstack<-raster::raster(paste0(dir_bm,"/",sp,"/",p,"/individual_projections/",p,"_",sp,"_ensemble.grd"))
      
    }
    
  }
}


# get ensemble w mean ensemble projections and w mean ensemble binaries
for (sp in sp.n){
  for (p in projects){
    
    # average weighted mean ensemble projections
    emstack<-raster::stack(paste0(dir_bm,"/",sp,"/",p,"/",p,"_",sp,"_ensemble.grd"))
    emwmeans <- raster::subset(emstack, grep('EMwmean', names(emstack), value = T))
    emwmean<-raster::stackApply(emwmeans,indices=c(rep(1,raster::nlayers(emwmeans))),fun=mean)
    emwmean<-emwmean/10
    raster::plot(emwmean,main=paste0(sp," Weighted Mean Ensemble for ",p),xlab="Longitude",ylab="Latitude")
    dir.create(paste0(dir_bm,"/weighted_mean_ensembles/"))
    raster::writeRaster(emwmean,file=paste0(dir_bm,"/weighted_mean_ensembles/",sp.n,"_",p,"_em-wmean.grd"),format="raster",overwrite=T)
    
    # get binary of average of w mean ensemble grds and write to file
    threshold<-"70"
    bin<-biomod2::BinaryTransformation(emwmean,threshold)
    raster::plot(bin)
    dir.create(paste0(dir_bm,"/binaries"))
    raster::writeRaster(bin,file=paste0(dir_bm,"/binaries/",sp.n,"_",p,"_",threshold,"_thresh_bin-wmean.grd"),format="raster",overwrite=T)
    
  }
}

# get range size change
projects
for (sp in sp.n){
  rangechange<-df()
  currentPred <- raster::raster(paste0(dir_bm,"/binaries/",sp.n,"_proj_current_",threshold,"_thresh_bin-wmean.grd"))
  f50Pred <- raster::raster(paste0(dir_bm,"/binaries/",sp.n,"_proj_rcp85_50_",threshold,"_thresh_bin-wmean.grd"))
  f70Pred<- raster::raster(paste0(dir_bm,"/binaries/",sp.n,"_proj_rcp85_70_",threshold,"_thresh_bin-wmean.grd"))
  
  sp::spplot(currentPred)
  sp::spplot(f50Pred)
  sp::spplot(f70Pred)
  # call the Range size function for 2050 average
  rangechange50 <- biomod2::BIOMOD_RangeSize(
    CurrentPred=currentPred,
    FutureProj=f50Pred)
  
  # see the results
  rangechange50$Compt.By.Models
  df50chg<-as.data.frame(rangechange70$Compt.By.Models)
  raster::plot(rangechange50$Diff.By.Pixel)
  rangechange70 <- biomod2::BIOMOD_RangeSize(
    CurrentPred=currentPred,
    FutureProj=f70Pred)
  
  # same for 2070 average
  rangechange70$Compt.By.Models
  df70chg<-as.data.frame(rangechange70$Compt.By.Models)
  df70chg
  raster::plot(rangechange70$Diff.By.Pixel)
  rangechange70
  
}


#################################################################
# CAPTURE PLOTS OF PROJECTIONS BY MODEL
#################################################################

# write individual model to image grouped by model
dir_projs<-paste0(dir_figs,"/",myRespName,"/projs")
dir.create(dir_projs,recursive=T)

print(paste0("Saving Individual Model Predictions by Model onto Current Data for ",myRespName))


eval_metrics<-c( 'KAPPA', 'TSS')
allmodels<-c("GLM","GAM","GBM","ANN","CTA","RF","MARS","FDA",'MAXENT.Phillips','MAXENT.Tsuruoka')

mod_proj <- get_predictions(myBiomodProj) 
for (mod in allmodels){
  indices<-grep(paste0("*",mod), layerNames(mod_proj))
  modprojs<-subset(mod_proj, indices, drop = TRUE)
  png(filename=paste0(dir_projs,"/",myRespName,"_",mod,"_projections.png"))
  plot(modprojs)
  dev.off()
}

print(paste0("Saving Individual Model Predictions by Model onto Future (2070) Data for ",myRespName))

f_70_proj <- get_predictions(myBiomodProjFuture70) 
for (mod in allmodels){
  indices<-grep(paste0("*",mod), layerNames(mod_proj))
  modprojs<-subset(f_70_proj, indices, drop = TRUE)
  png(filename=paste0(dir_projs,"/",myRespName,"_",mod,"_rcp85_70_projections.png"))
  plot(modprojs)
  dev.off()
}

print(paste0("Saving Individual Model Predictions by Model onto Future (2050) Data for ",myRespName))

f_50_proj <- get_predictions(myBiomodProjFuture50) 
for (mod in allmodels){
  indices<-grep(paste0("*",mod), layerNames(mod_proj))
  modprojs<-subset(f_50_proj, indices, drop = TRUE)
  png(filename=paste0(dir_projs,"/",myRespName,"_",mod,"_rcp85_50_projections.png"))
  plot(modprojs)
  dev.off()
}




# plot ensemble forecasts grouped by metric
projects<-c("proj_current","proj_rcp85_50","proj_rcp85_70")
eval_metrics<-c( 'KAPPA', 'TSS')


current_em_proj <- get_predictions(myBiomodEF) 

print(paste0("Saving Current Ensemble Plots for  ",myRespName))

for (eval in eval_metrics){
  dir_forecasts<-paste0(dir_figs,"/",myRespName,"/forecasts")
  dir.create(dir_forecasts,recursive=T)
  indices<-grep(paste0("*",mod), layerNames(mod_proj))
  forecasts<-subset(current_em_proj, indices, drop = TRUE)
  png(filename=paste0(dir_forecasts,"/",myRespName,"_",eval,"_em_forecasts.png"))
  plot(forecasts)
  dev.off()
}

print(paste0("Saving Future (2070) Ensemble Plots for  ",myRespName))

f70forecasts <- get_predictions(f70BiomodEF) 

for (eval in eval_metrics){
  dir_forecasts<-paste0(dir_figs,"/",myRespName,"/forecasts")
  dir.create(dir_forecasts,recursive=T)
  indices<-grep(paste0("*",mod), layerNames(mod_proj))
  forecasts<-subset(f70forecasts, indices, drop = TRUE)
  png(filename=paste0(dir_forecasts,"/",myRespName,"_",eval,"_em_f70_forecasts.png"))
  plot(forecasts)
  dev.off()
}

print(paste0("Saving Future (2050) Ensemble Plots for  ",myRespName))

f50forecasts <- get_predictions(f50BiomodEF) 

for (eval in eval_metrics){
  dir_forecasts<-paste0(dir_figs,"/",myRespName,"/forecasts")
  dir.create(dir_forecasts,recursive=T)
  indices<-grep(paste0("*",mod), layerNames(mod_proj))
  forecasts<-subset(f50forecasts, indices, drop = TRUE)
  png(filename=paste0(dir_forecasts,"/",myRespName,"_",eval,"_em_f50_forecasts.png"))
  plot(forecasts)
  dev.off()
}
# myBiomodEF
# EF_70stack<- raster::stack(paste0(dir_bm,"/",sp.n,"/proj_current/proj_rcp85_70",sp.n,"_ensemble.grd"))
# plot a layer in the ensemble stack
# plot(EF_stack[[6]])

# myBiomodEF
# EF_50stack<- raster::stack(paste0(dir_bm,"/",sp.n,"/proj_current/proj_rcp85_50",sp.n,"_ensemble.grd"))
# plot a layer in the ensemble stack
# plot(EF_stack[[6]])


#################################################################
# CAPTURE MODEL RESPONSE CURVES TO FILE
#################################################################
print(paste0("Writing Response Curve Plots to File for ",myRespName))

myBiomodModelOut
mod<-models[4]
var <- sp.n[1]
for (var in sp.n){
  
  for (mod in models){
    modelfile<-load(paste0(dir_bm,"/",var,"/",var,".",var,"_current.models.out"))
    modelsubset<-get(modelfile)
    modelsubset
    response.plot2(models= loadedmodel,
                   Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                   show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
                   save.file = "tiff",
                   name=paste0(dir_curves,"/",myRespName,"_",mod,"_curves"),
                   col = c("blue", "red"),
                   legend = TRUE,
                   data_species = get_formal_data(myBiomodModelOut,'resp.var'),
                   ImageSize=1000)
    
  }
}

for (mod in models){
  for (sp in sp.n){
    
    dir_curves<-paste0(dir_figs,"/",sp,"/response-curves")
    dir.create(dir_curves,recursive=T)
    dir.create(paste0(dir_figs,"/",sp))
    modelfile<-load(paste0(dir_bm,"/",sp,"/",sp,".",sp,"_current.models.out"))
    myBiomodModelOut<-get(modelfile)
    loadedmodel<-biomod2::BIOMOD_LoadModels(myBiomodModelOut,models=mod)
    dir.create(paste0(dir_figs,"/",myRespName))
    
    
    response.plot2(models= loadedmodel,
                   Data = get_formal_data(myBiomodModelOut,'expl.var'),
                   show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
                   save.file = "tiff",
                   name=paste0(dir_curves,"/",myRespName,"_",mod,"_curves"),
                   col = c("blue", "red"),
                   legend = TRUE,
                   data_species = get_formal_data(myBiomodModelOut,'resp.var'))
  }
  
  #################################################################
  # GET MODEL SCORES GRAPH
  #################################################################
  print(paste0("Capturing Model Scores Graph for ",myRespName))
  
  dir.create(paste0(dir_figs,"/",myRespName))
  png(filename=paste0(dir_figs,"/",myRespName,"/",myRespName,"_model_scores-kappa-tss.png"))
  models_scores_graph(myBiomodModelOut,metrics = c( 'KAPPA', 'TSS'),by = 'models',plot = TRUE)
  dev.off()
  
  #################################################################
  # GET PLOT OF OBSERVATIONS
  #################################################################
  print(paste0("Plotting Observations for ",myRespName))
  
  png(filename=paste0(dir_figs,"/",myRespName,"/",myRespName,"_observations.png"))
  presab<-myResp
  presab[is.na(presab)]<-0
  level.plot(data.in = presab,
             XY = myRespXY,
             color.gradient = "red",
             cex = .7,
             level.range = c(min(presab), max(presab)),
             show.scale = FALSE,
             title = paste0(myRespName," Observations"),
             SRC=FALSE,
             ImageSize="large"
  )
  
  dev.off()
  # for (p in projs){
  #  ensemblestack<-raster::stack(paste0(dir_bm,"/",sp.n,"/",p,"/",p,"_",sp.n,"_ensemble.grd"))
  #  ensembleconsensus<-stackApply(ensemblestack,indices=rep(1:1,nlayers(ensemblestack)),fun=mean)
  #  ensembleconsensus100<-ensembleconsensus/1000
  #  writeRaster(ensembleconsensus100,file=paste0(dir_bm,"/",sp.n,"/",p,"/",p,"_",sp.n,"_mean_consensus.grd"),format="raster",overwrite=T)
  #}
  
  # define a mask of studied
  alphaMap <- reclassify(subset(myExpl,1), c(-Inf,Inf,0))
  currentensemble<-stack(file.path(variety,
                                   "proj_current", 
                                   paste("proj_current_",
                                         variety, 
                                         "_ensemble.grd", sep="")))
  names(currentensemble)
  
  
  
  
  # get alpha diversity map
  pa<-read.csv(file=paste0(dir_out,"/pa_dataframe.csv"))
  pa<-data.frame(pa)
  
  # get vector of species names
  names<-paste0(colnames(pa))
  varieties <- dput(names[c(49)])
  
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
  