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

#library(mgcv)


#################################################################
# PRELIMINARY DATA FORMATTING
#################################################################



# read in model raster stacks in LONG/LAT
# can't project onto equal area because Maxent needs resolution with equal x, y
# this data came from lat/long, so the projected 
library(raster)
library(rgdal)
library(biomod2)


#################################################################
# INITIALIZE FUNCTION TO APPLY TO EACH VARIETY
#################################################################



options(max.print=1000000)  # set max.print option high to capture outputs


setwd(dir_bm) 

# get species vector name by getting all with ensemble projections for 2070 time period
# or whichever was last in the BioModApply function/loop

emprojfiles<-Sys.glob(paste0(dir_bm,"/*/*/","*70.ensemble.projection.out")) # get ensemble files out
emproj.dirs<-sub(paste0(dir_bm,"/"),"",emprojfiles) # remove everything from working directory path
sp.n<-sub(" */.*", "",emproj.dirs) # remove everything before first slash to get variety names that have ensemble models
sp.n

#################################################################
# get variable importance csvs and merge into one
#################################################################

# under construction 

# following https://www.r-bloggers.com/merging-multiple-data-files-into-one-data-frame/
#var_imp_files<-Sys.glob(paste0(dir_out,"/*/*_var_imp.csv")) # get ensemble files out

#var_imp.df = lapply(var_imp_files, function(x){read.csv(file=x,header=T)})
#View(var_imp.df)
setwd(dir_bm) 

#library(data.table)
#multmerge = function(filenames){
# rbindlist(lapply(filenames, function(x){read.csv(x, stringsAsFactors = F, sep=',',header = TRUE)

#}))

#}


#################################################################
# CAPTURE MODEL RESPONSE CURVES TO FILE
#################################################################

setwd(dir_bm) 
#options(max.print=1000000)  # set max.print option high to capture outputs
#path.to.maxent.jar<-file.path(getwd(),"maxent.jar") # define maxent jar location



models = c("MAXENT.Phillips")
metrics = c( 'KAPPA', 'TSS', 'ROC')
projects<-c("proj_current","proj_rcp85_50","proj_rcp85_70")

models = c("MAXENT.Phillips")
mod<-models

for (sp.n in sp.n){
  #PostBioApply <-function(sp.n) {
  #tryCatch({

  setwd(dir_bm)

  #for (mod in models){

  #################################################################
  # CAPTURE MODEL RESPONSE CURVES TO FILE BY ALGORITHM
  #################################################################




  #################################################################
  # GET MODEL EVALUATIONS
  #################################################################


  modelfile<-load(paste0(dir_bm,"/",sp.n,"/",sp.n,".",sp.n,"_current.models.out"))
  myBiomodModelOut<-get(modelfile)
  myBiomodModelOut

  #print(paste0("Capturing Model Scores Graph for ",sp.n))
  #dir_evals<-paste0(dir_out,"/model-evals/")
  #dir.create(dir_evals,recursive=TRUE,showWarnings=FALSE)


  # write em models built
  #capture.output(get_built_models (myBiomodModelOut),
  #               file=paste0(dir_evals,"/",sp.n,"_models.txt"))

  # capture em model evals
  #print(paste0("Capturing Models Evaluations ",sp.n))


  #evalmods<-get_evaluations(myBiomodModelOut,as.data.frame=TRUE)
  #evalmods$variety<-sp.n
  #write.csv(evalmods,file=paste0(dir_evals,"/",sp.n,"_models_eval.csv"))

  #do.call(file.remove,list(list.files(pattern="temp*")))

  ### eval current model

  #print(paste0("Capturing Model Ensemble Evaluations for ",sp.n))
  #evalmods<-get_evaluations(myBiomodModelOut,as.data.frame=TRUE)
  #write.csv(evalmods,file=paste0(dir_evals,"/",sp.n,"_evals-df.csv"))

  #################################################################
  # GET ENSEMBLE MODEL EVALUATIONS
  #################################################################

  emmodelfile<-load(paste0(dir_bm,"/",sp.n,"/",sp.n,".",sp.n,"_currentensemble.models.out"))
  myBiomodEM<-get(emmodelfile)
  myBiomodEM
  #print(paste0("Capturing Model Scores Graph for ",sp.n))
  #dir_enevals<-paste0(dir_out,"/ensemble-evals/")
  #dir.create(dir_enevals,recursive=TRUE,showWarnings=FALSE)


  #print(paste0("Capturing Ensemble Model Outputs for  ",sp.n))

  # write em models built
  #capture.output(get_built_models (myBiomodEM),
  #               file=paste0(dir_enevals,"/",sp.n,"_em_models.txt"))

  # capture em model evals
  #print(paste0("Capturing Ensemble Models Evaluations ",sp.n))
  #evalmodsEM<-get_evaluations(myBiomodEM,as.data.frame=TRUE)
  #evalmodsEM$variety<-sp.n
  #write.csv(evalmodsEM,file=paste0(dir_enevals,"/",sp.n,"_models_EM_eval.csv"))

  #do.call(file.remove,list(list.files(pattern="temp*")))

  ### eval current model

  #print(paste0("Capturing Model Ensemble Evaluations for ",sp.n))
  #enevalmods<-get_evaluations(myBiomodEM,as.data.frame=TRUE)
  #write.csv(enevalmods,file=paste0(dir_enevals,"/",sp.n,"_em_evals-df.csv"))
  
  #################################################################
  # GET VARIABLE IMPORTANCE
  #################################################################
  
  
  print(paste0("Capturing Model Scores Graph for ",sp.n))
  dir_var<-paste0(dir_out,"/var-importance/")
  
  dir.create(dir_var,recursive=TRUE,showWarnings=FALSE)
  
  
  modvarimp<-get_variables_importance(myBiomodModelOut)
  modvarimp.df<-as.data.frame(modvarimp)
  modvarimp.df$variety<-sp.n
  write.csv(modvarimp.df,file=paste0(dir_var,"/",sp.n,"_model-var-imp-df.csv"))
  
  EMvarimp<-get_variables_importance(myBiomodEM)
  EMvarimp.df<-as.data.frame(rowMeans(EMvarimp))
  EMvarimp.df$variety<-sp.n
  
  write.csv(EMvarimp.df,file=paste0(dir_var,"/",sp.n,"_EM-mean-var-imp-df.csv"))
  
  #################################################################
  # CAPTURE MODEL RESPONSE CURVES TO FILE BY MODEL ALGORITHM RUN
  #################################################################
# 
#   print(paste0("Capturing Model Scores Graph for ",sp.n))
#   dir_curves<-paste0(dir_out,"/response-curves/")
#   dir.create(dir_curves,recursive=TRUE,showWarnings=FALSE)
# 
#   modelsubset<-BIOMOD_LoadModels(myBiomodModelOut,models=mod)
# 
#   pdf(width = 6, height = 6)
#   response.plot2(models=modelsubset,
#                  Data = get_formal_data(myBiomodModelOut,'expl.var'),
#                  show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
#                  save.file = "pdf",
#                  col=rainbow(length(myBiomodModelOut@models.computed)),
#                  name=paste0(dir_curves,sp.n,"_",mod,"_curves"),
#                  do.bivariate=FALSE,
#                  fixed.var.metric = 'median',
#                  data_species = get_formal_data(myBiomodModelOut,'resp.var'),
#                  ImageSize=5000,
#                  plot=TRUE)


  ### 3D Response curves loop by individual model run
  #for (i in 1:length(modelsubset)){

  #myRespPlot3D <- response.plot2(models  = modelsubset[i],
  #Data = get_formal_data(myBiomodModelOut,'expl.var'),
  #show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
  #do.bivariate = TRUE,
  #fixed.var.metric = 'median',
  #data_species = get_formal_data(myBiomodModelOut,'resp.var'),
  #display_title=FALSE,
  #col = c("blue", "red"),
  #legend = TRUE,
  #save.file = "tiff",
  #name=paste0(dir_curves,sp.n,"_",mod,"_",i,"_3d_curves"),
  #ImageSize=500,
  # plot=TRUE)


  # } # response plot3d loop by individual model



  #################################################################
  # GET MODEL SCORES GRAPH   -- under construction
  #################################################################



  #print(paste0("Capturing Model Scores Graph for ",sp.n))
  #dir_modscores<-paste0(dir_out,"/model-score-graphs/")
  #dir.create(dir_modscores,recursive=TRUE,showWarnings=FALSE)


  # for (metrics in seq_along(metrics)) ....

  #pdf(filename=paste0(dir_modscores,"/",sp.n,"_model_scores-",paste(metrics[1:2], collapse = "-"),".pdf"))
  #models_scores_graph(myBiomodModelOut,metrics = metrics[1:2],plot = TRUE)
  #dev.off()


  #pdf(filename=paste0(dir_modscores,"/",sp.n,"_model_scores-",paste(metrics[2:3], collapse = "-"),".pdf"))
  #models_scores_graph(myBiomodModelOut,metrics = metrics[2:3],by = 'cv_run',plot = TRUE, width = 12, height = 17, family = "Helvetica")
  #dev.off()


  #pdf(filename=paste0(dir_modscores,"/",sp.n,"_model_scores-",paste(metrics[c(1,3)], collapse = "-"),".pdf"))
  #models_scores_graph(myBiomodModelOut,metrics = metrics[c(1,3)],by = 'cv_run',plot = TRUE, width = 12, height = 17, family = "Helvetica")
  #dev.off()



  # get ensemble w mean ensemble projections and w mean ensemble binaries

  #################################################################
  # GET WEIGHTED MEAN BINARIES WITH THERSHOLD
  #################################################################
# 
#   for (p in projects){
# 
#     dir_wmean<-paste0(dir_out,"/wmean-consensus")
#     dir.create(dir_wmean,recursive=TRUE,showWarnings=FALSE)
# 
#     # average weighted mean ensemble projections
#     emstack<-raster::stack(paste0(dir_bm,"/",sp.n,"/",p,"/",p,"_",sp.n,"_ensemble.grd"))
#     emwmeans <- raster::subset(emstack, grep('EMwmean', names(emstack), value = T))
# 
#     emwmean<-raster::stackApply(emwmeans,indices=c(rep(1,raster::nlayers(emwmeans))),fun=mean)
#     #raster::plot(emwmean,main=paste0(sp.n," Weighted Mean Ensemble for ",p),xlab="X",ylab="Y")
#     raster::writeRaster(emwmean,file=paste0(dir_wmean,"/",sp.n,"_",p,"_em-wmean.grd"),format="raster",overwrite=T)
# 
# 
#     # get binary of average of w mean ensemble grds and write to file
# 
#     rc <- function(x) {
# 
#       ifelse(x <=  750, 0,
# 
#              ifelse(x >  0.75, 1, NA)) }
# 
#     bin <- calc(emwmean, fun=rc)
#     threshold<-"75"
#     # bin<-biomod2::BinaryTransformation(emwmean,threshold)
#     #raster::plot(bin)
# 
#     dir_binary<-paste0(dir_out,"/binaries")
#     dir.create(dir_binary,recursive=TRUE,showWarnings=FALSE)
# 
#     avg_wmean_bin_grd<-paste0(dir_binary,"/",sp.n,"_",p,"_",threshold,"_thresh_bin-avg-wmean.grd")
#     raster::writeRaster(bin,file=avg_wmean_bin_grd,format="raster",overwrite=T)
# 
#   }# for p in projects

 #################################################################
  # USE BINARY TO CACULATE RANGE SIZE CHANGE 
  #################################################################
  # dir_binary<-paste0(dir_out,"/binaries")
  # 
  # dir_range<-paste0(dir_out,"/range-size")
  # 
  # dir_rangedf<-paste0(dir_out,"/range-size/dfs")
  # 
  # dir_rangeras<-paste0(dir_out,"/range-size/rasters")
  # 
  # dir_rangeR<-paste0(dir_out,"/range-size/RData")
  #     threshold<-"75"
  # 
  # currentPred <- raster::raster(paste0(dir_binary,"/",sp.n,"_proj_current_",threshold,"_thresh_bin-avg-wmean.grd"))
  # f50Pred <-  raster::raster(paste0(dir_binary,"/",sp.n,"_proj_rcp85_50_",threshold,"_thresh_bin-avg-wmean.grd"))
  # f70Pred<-  raster::raster(paste0(dir_binary,"/",sp.n,"_proj_rcp85_70_",threshold,"_thresh_bin-avg-wmean.grd"))
  # 
  # f50Pred<-resample(f50Pred, currentPred, method="ngb")
  # f70Pred<-resample(f70Pred, currentPred, method="ngb")
  # 
  # #sp::spplot(currentPred)
  # #sp::spplot(f50Pred)
  # #sp::spplot(f70Pred)
  # 
  # #################################################################
  # # CHANGE BY 2050
  # #################################################################
  # 
  # # call the Range size function for 2050 average
  # rangechange50 <- biomod2::BIOMOD_RangeSize(
  #   CurrentPred=currentPred,
  #   FutureProj=f50Pred)
  # 
  # # see the results
  # rangechange50$Compt.By.Models
  # df50chg<-as.data.frame(rangechange50$Compt.By.Models)
  # df50chg$time<-"2041-2060"
  # df50chg$species<-sp.n
  # df50chg
  # 
  # write.csv(df50chg,file=paste0(dir_rangedf,"/",sp.n,"_f50_rangesize.csv"))
  # writeRaster(rangechange50$Diff.By.Pixel ,file=paste0(dir_rangeras,"/",sp.n,"_2050_diffbypixel.grd"),overwrite=TRUE)    
  # save(rangechange50,file=paste0(dir_rangeR,"/",sp.n,"_2050_rangechange.Rdata"))
  # 
  # 
  # #################################################################
  # # CHANGE BY 2070
  # #################################################################
  # 
  # # call the Range size function for 2050 average
  # rangechange70 <- biomod2::BIOMOD_RangeSize(
  #   CurrentPred=currentPred,
  #   FutureProj=f70Pred)
  # 
  # # see the results
  # rangechange70$Compt.By.Models
  # df70chg<-as.data.frame(rangechange70$Compt.By.Models)
  # df70chg$time<-"2061-2080"
  # df70chg$species<-sp.n
  # df70chg
  # 
  # write.csv(df50chg,file=paste0(dir_rangedf,"/",sp.n,"_f70_rangesize.csv"))
  # writeRaster(rangechange70$Diff.By.Pixel ,file=paste0(dir_rangeras,"/",sp.n,"_2070_diffbypixel.grd"),overwrite=TRUE)    
  # save(rangechange70,file=paste0(dir_rangeR,"/",sp.n,"_2070_rangechange.Rdata"))
  # 
  #}
}

# threshold<-"75"
# 
# for (sp.n in sp.n){
#  
#   
#   
# } # PostBioApply function 
# 
# 



#################################################################
# GET ALPHA DIVERSITY
#################################################################

# # following https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/biomod2/inst/doc/Multi_species_computation.pdf?revision=315&root=biomod


# for (p in projects){
#   dir_range<-paste0(dir_out,"/range-size")
#   bin_grds<-Sys.glob(paste0(dir_binary,"/*",p,"*.grd")) # get ensemble files out
#   bin_stack<-stack(bin_grds)
#   alpha_div<-raster::stackApply(bin_stack,indices=c(rep(1,raster::nlayers(bin_stack))),fun=sum)
#   writeRaster(alpha_div ,file=paste0(dir_range,"/",p,"_alpha_diversity.grd"),overwrite=TRUE)  
# }






























#################################################################
# CAPTURE PLOTS OF PROJECTIONS BY MODEL
#################################################################

# write individual model to image grouped by model
#dir_projs<-paste0(dir_figs,"/",myRespName,"/projs")
#dir.create(dir_projs,recursive=T)

#print(paste0("Saving Individual Model Predictions by Model onto Current Data for ",myRespName))


#mod_proj <- get_predictions(myBiomodProj) 
#for (mod in allmodels){
#  indices<-grep(paste0("*",mod), layerNames(mod_proj))
#  modprojs<-subset(mod_proj, indices, drop = TRUE)
#  tiff(filename=paste0(dir_projs,"/",myRespName,"_",mod,"_projections.tiff", height = 12, width = 17, units = 'cm', 
#                       compression = "lzw", res = 300))
#  plot(modprojs)
#  dev.off()
#}

#print(paste0("Saving Individual Model Predictions by Model onto Future (2070) Data for ",myRespName))

#f_70_proj <- get_predictions(myBiomodProjFuture70) 
#for (mod in allmodels){
#  indices<-grep(paste0("*",mod), layerNames(mod_proj))
#  modprojs<-subset(f_70_proj, indices, drop = TRUE)
#  tiff(filename=paste0(dir_projs,"/",myRespName,"_",mod,"_rcp85_70_projections.tiff", width = 12, height = 17, family = "Helvetica"))
#  plot(modprojs)
#  dev.off()
#}

#print(paste0("Saving Individual Model Predictions by Model onto Future (2050) Data for ",myRespName))

#f_50_proj <- get_predictions(myBiomodProjFuture50) 
#for (mod in allmodels){
#  indices<-grep(paste0("*",mod), layerNames(mod_proj))
#  modprojs<-subset(f_50_proj, indices, drop = TRUE)
#  tiff(filename=paste0(dir_projs,"/",myRespName,"_",mod,"_rcp85_50_projections.tiff", width = 12, height = 17, family = "Helvetica"))
#  plot(modprojs)
#  dev.off()
#}



#################################################################
# merge variable importance csvs
#################################################################

setwd(paste0(dir_out,"/var-importance/"))
emvarimp.files <- list.files(pattern="*EM-mean*")
emvarimp.paths<-paste0(getwd(),"/",emvarimp.files)
var.imp.merge<-data.frame()


### merge variable importance (ensemble) csvs into dataframe
merged.em.varimp <- 
  do.call(rbind,
          lapply(emvarimp.paths, read.csv))

### group variety and variable
dply<-merged.em.varimp %>% group_by_("X","variety")

library(tidyr)

### spread to get columns for each variety
tidy<-dply %>% spread(variety,rowMeans.EMvarimp.)

### transpose to get variables as columns
tidy.t<-t(tidy)
tidy.df<-as.data.frame(tidy.t)

### get variables as column names (from first row)
colnames(tidy.df) <- as.character(unlist(tidy.df[1,]))
tidy.df = tidy.df[-1, ]
View(tidy.df)



#################################################################
# merge model evaluation csvs
#################################################################

setwd(paste0(dir_out,"/ensemble-evals/"))
emval.files <- list.files(pattern="*EM_eval*")
emval.paths<-paste0(getwd(),"/",emval.files)

### merge variable importance (ensemble) csvs into dataframe
merged.em.evals <- 
  do.call(rbind,
          lapply(emval.paths, read.csv))

View(merged.em.evals)
colnames(merged.em.evals)


### group variety and variable
dply<-merged.em.evals %>% group_by(variety,Model.name,Eval.metric) %>% summarize(Testing.data=mean(Testing.data),Sensitivity=mean(Sensitivity),Specificity=mean(Specificity)) 
View(dply)

var.avg<-merged.em.evals %>% group_by(variety,Eval.metric) %>% summarize(Testing.data=mean(Testing.data),Sensitivity=mean(Sensitivity),Specificity=mean(Specificity)) 
View(var.avg)

spread(dply,Eval.metric,Testing.data)
### spread to get columns for each variety
tidy<-dply %>% spread(Eval.metric,variety)
View(tidy)


#################################################################
# diversity zonal stats
#################################################################
#import required libraries
library(maptools)
library(raster)
library(rgdal)
par(mar = rep(2, 4))

dir_range<-paste0(dir_out,"/range-size")
#list files (in this case raster TIFFs)
alpha.div <- raster(paste0(dir_range,"/proj_current_alpha_diversity.grd"))


#read-in the polygon shapefile
poly <- readOGR("E:\\thesis\\01_data\\ind\\presindigw\\presindigw.shp",layer="presindigw")
ob <- SpatialPolygons(poly@polygons,proj4string=poly@proj4string)
ob<-spTransform(ob,proj4string(alpha.div))
poly<-spTransform(poly,proj4string(alpha.div))

#extract raster cell count (sum) within each polygon area (poly)
ex <- raster::extract(alpha.div,poly,fun=median,na.rm=TRUE,df=TRUE)
poly@data$pobpct<-poly@data$pobindi/poly@data$pobtot
df<-merge(ex, poly@data)

library(corrgram)


cor.pearson <-cor.test(df$index_1,df$pobpct,  method = "pearson")
cor.pearson
cor.kendall <-cor.test(df$index_1,df$pobpct,  method = "kendall ")
cor.kendall
cor.spearman <-cor.test(df$index_1,df$pobpct,  method = "spearman ")
cor.spearman

#write to a CSV file
write.csv(df, file = "./path/to/data/CSV.csv")
