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

# 
#   modelfile<-load(paste0(dir_bm,"/",sp.n,"/",sp.n,".",sp.n,"_current.models.out"))
#   myBiomodModelOut<-get(modelfile)
#   myBiomodModelOut

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
# 
#   emmodelfile<-load(paste0(dir_bm,"/",sp.n,"/",sp.n,".",sp.n,"_currentensemble.models.out"))
#   myBiomodEM<-get(emmodelfile)
#   myBiomodEM
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
  # 
  # 
  # print(paste0("Capturing Model Scores Graph for ",sp.n))
  # dir_var<-paste0(dir_out,"/var-importance/")
  # 
  # dir.create(dir_var,recursive=TRUE,showWarnings=FALSE)
  # 
  # 
  # modvarimp<-get_variables_importance(myBiomodModelOut)
  # modvarimp.df<-as.data.frame(modvarimp)
  # modvarimp.df$variety<-sp.n
  # write.csv(modvarimp.df,file=paste0(dir_var,"/",sp.n,"_model-var-imp-df.csv"))
  # 
  # EMvarimp<-get_variables_importance(myBiomodEM)
  # EMvarimp.df<-as.data.frame(rowMeans(EMvarimp))
  # EMvarimp.df$variety<-sp.n
  # 
  # write.csv(EMvarimp.df,file=paste0(dir_var,"/",sp.n,"_EM-mean-var-imp-df.csv"))
  # 
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
  dir_binary<-paste0(dir_out,"/binaries")

  dir_range<-paste0(dir_out,"/range-size")

  dir_rangedf<-paste0(dir_out,"/range-size/dfs")

  dir_rangeras<-paste0(dir_out,"/range-size/rasters")

  dir_rangeR<-paste0(dir_out,"/range-size/RData")
      threshold<-"75"

  currentPred <- raster::raster(paste0(dir_binary,"/",sp.n,"_proj_current_",threshold,"_thresh_bin-avg-wmean.grd"))
  f50Pred <-  raster::raster(paste0(dir_binary,"/",sp.n,"_proj_rcp85_50_",threshold,"_thresh_bin-avg-wmean.grd"))
  f70Pred<-  raster::raster(paste0(dir_binary,"/",sp.n,"_proj_rcp85_70_",threshold,"_thresh_bin-avg-wmean.grd"))

  f50Pred<-resample(f50Pred, currentPred, method="ngb")
  f70Pred<-resample(f70Pred, currentPred, method="ngb")

  #sp::spplot(currentPred)
  #sp::spplot(f50Pred)
  #sp::spplot(f70Pred)

  #################################################################
  # CHANGE BY 2050
  #################################################################

  # call the Range size function for 2050 average
  rangechange50 <- biomod2::BIOMOD_RangeSize(
    CurrentPred=currentPred,
    FutureProj=f50Pred)

  # see the results
  rangechange50$Compt.By.Models
  df50chg<-as.data.frame(rangechange50$Compt.By.Models)
  df50chg$time<-"2041-2060"
  df50chg$species<-sp.n
  df50chg

  write.csv(df50chg,file=paste0(dir_rangedf,"/",sp.n,"_f50_rangesize.csv"))
  writeRaster(rangechange50$Diff.By.Pixel ,file=paste0(dir_rangeras,"/",sp.n,"_2050_diffbypixel.grd"),overwrite=TRUE)
  save(rangechange50,file=paste0(dir_rangeR,"/",sp.n,"_2050_rangechange.Rdata"))


  #################################################################
  # CHANGE BY 2070
  #################################################################

  # call the Range size function for 2050 average
  rangechange70 <- biomod2::BIOMOD_RangeSize(
    CurrentPred=currentPred,
    FutureProj=f70Pred)

  # see the results
  rangechange70$Compt.By.Models
  df70chg<-as.data.frame(rangechange70$Compt.By.Models)
  df70chg$time<-"2061-2080"
  df70chg$species<-sp.n
  df70chg

  write.csv(df70chg,file=paste0(dir_rangedf,"/",sp.n,"_f70_rangesize.csv"))
  writeRaster(rangechange70$Diff.By.Pixel ,file=paste0(dir_rangeras,"/",sp.n,"_2070_diffbypixel.grd"),overwrite=TRUE)
  save(rangechange70,file=paste0(dir_rangeR,"/",sp.n,"_2070_rangechange.Rdata"))

  }
#}

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

# dir_binary<-paste0(dir_out,"/binaries")
# 
#  for (p in projects){
#    dir_range<-paste0(dir_out,"/range-size")
# 
# dir_binary<-paste0(dir_out,"/binaries")
# 
#    bin_grds<-Sys.glob(paste0(dir_binary,"/*",p,"*.grd")) # get ensemble files out
#    bin_stack<-stack(bin_grds)
#    alpha_div<-raster::stackApply(bin_stack,indices=c(rep(1,raster::nlayers(bin_stack))),fun=sum)
#    writeRaster(alpha_div ,file=paste0(dir_range,"/",p,"_alpha_diversity.grd"),overwrite=TRUE)  
#  }
# 
# 
# library(rasterVis)
# library(maps)
# library(mapdata)
# library(maptools)
# library(sp)
# 
dir_wmean<-paste0(dir_out,"/wmean-consensus")
dir.create(dir_wmean,recursive=TRUE,showWarnings=FALSE)
dir_binary<-paste0(dir_out,"/binaries")
dir.create(dir_binary,recursive=TRUE,showWarnings=FALSE)
dir_rangeras<-paste0(dir_out,"/range-size/rasters")
dir_range<-paste0(dir_out,"/range-size")
dir_rangedf<-paste0(dir_out,"/range-size/dfs")
dir_rangeR<-paste0(dir_out,"/range-size/RData")
sp<-sp.n[1]
    for (sp in sp.n){
      gc()
        dir_wmean<-paste0(dir_out,"/wmean-consensus")

        wmean<-Sys.glob(paste0(dir_wmean,"/",sp,"*.grd")) # get ensemble files out
        current.wmean<-raster(wmean[1])
        f50wmean<-resample(raster(wmean[2]),current.wmean)
        f70wmean<-resample(raster(wmean[2]),current.wmean)
        wmstack<-stack(current.wmean,f50wmean,f70wmean)
        names(wmstack)<-c("1970-2000","2041-2060","2061-2080")
        mapTheme <- rasterTheme(region=brewer.pal(8,"Greens"))

        plt <-levelplot(wmstack, margin = list(FUN = 'median'), par.settings=mapTheme,
        xlab="X", ylab="Y",main=paste0(sp," \nProjected Relative Probability of Presence"), names.attr=c("1970-2000","2041-2060","2061-2080"))

        pdf(paste0(dir_wmean,"/",sp,"-wmean-em-projs.pdf"))
        print(plt  )
        dev.off()

        bin_grds<-Sys.glob(paste0(dir_binary,"/",sp,"*.grd"))
        current<-raster(bin_grds[1])
        f50bin<-resample(raster(bin_grds[2]),current)
        f70bin<-resample(raster(bin_grds[2]),current)

        binstack<-stack(current,f50bin,f70bin)

        binplt <-levelplot(binstack, par.setqtings=mapTheme,xlab="X", ylab="Y",main=paste0(sp," Binary Predictions (.75 threshold) \nfrom Average Weighted Mean Ensemble"),maxpixels= 20000000, names.attr=c("1970-2000","2041-2060","2061-2080"))

        pdf(paste0(dir_binary,"/",sp,"-wmean-bins.pdf"))
        print(binplt )
        dev.off()
        gc()
    }

library(gtable)
library(grid)
library(ggplot2)



## plot binary differences by pixel for both future periods

for (sp in sp.n){
  gc()
  chg<-Sys.glob(paste0(dir_rangeras,"/",sp,"*.grd")) # get ensemble files out
  chgras<-raster(chg[1])
  val <- getValues(chgras)
  xy <- as.data.frame(xyFromCell(chgras,1:ncell(chgras)))
  xy <- cbind(xy,val)
  
  cols<-c("-2"="firebrick", "-1"="skyblue2","0"="cornsilk","1"="royalblue4")
  p <-  ggplot(na.omit(xy), aes(x=x, y=y, fill=factor(val))) + 
    geom_raster()+ ggtitle(paste0(sp,"\nBinary Classification Change by Pixel"), subtitle ="2041 - 2060")+
    coord_equal() + scale_fill_manual(values = cols, name="",drop=FALSE,
                                      breaks = c("-2", "-1", "0","1"),labels=c("Loss", "Stable", "Neither","Migrate")) 
  
    change<-p 
  chgras<-raster(chg[2])
  val <- getValues(chgras)
  xy <- as.data.frame(xyFromCell(chgras,1:ncell(chgras)))
  xy <- cbind(xy,val)
  p <- ggplot(na.omit(xy), aes(x=x, y=y, fill=factor(val))) + 
    geom_raster() +ggtitle(paste0(""), subtitle ="2061 - 2080")+
    coord_equal() + scale_fill_manual(values = cols, drop=FALSE, name="",
                                      breaks = c("-2", "-1", "0","1"),labels=c("Loss", "Stable", "Neither","Migrate")) 
  
  change2<-p 
  library(gtable)
  library(grid)
  library(ggplot2)
  g2 <- ggplotGrob(change)
  g3 <- ggplotGrob(change2)
  g <- rbind(g2, g3, size = "first")
  g$widths <- unit.pmax(g2$widths, g3$widths)
  grid.newpage()
  
  ggsave(paste0(sp,"_change_by_pixels.pdf"), plot = grid.draw(g), device = "pdf", path = dir_rangeras,
         scale = 1)
}



setwd(dir_range)
alphadiv.files <- list.files(pattern="*.grd")
alphadiv.paths<-paste0(getwd(),"/",alphadiv.files)

div.stack<-stack(alphadiv.paths)
div.stack[[grep("current")]]div.stack[[2]]
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

View(merged.em.varimp)
library(dplyr)
### group variety and variable
dply<-merged.em.varimp %>% group_by_("X","variety")

library(tidyr)

### spread to get columns for each variety
tidy<-dply %>% spread(variety,rowMeans.EMvarimp.)

### transpose to get variables as columns
tidy.t<-t(tidy)
tidy.df<-as.data.frame(tidy.t)
colnames(tidy.df) <- as.character(unlist(tidy.df[1,]))
tidy.df= tidy.df[-1, ]
tidy.df<-tidy.df[,2:ncol(tidy.df)]
View(tidy.df)

write.csv(tidy.df,paste0(dir_out,"/var-importance/merged.csv"))

# Ward Hierarchical Clustering
colnames(tidy.df)
var.df<-tidy.df#[,c("IrrCult","Rain.fedCult","ind_pob_pcnt","Urban","Grass.Woodland","Water","vrm")]
var.df<-var.df[complete.cases(var.df),]
var.df<-var.df
d <- dist(var.df, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
rect.hclust(fit, k=5, border="red")
var.df$groups<-groups
colnames(var.df[,1:19])
for(i in c(1:19)) {
  var.df[,i] <- as.numeric(var.df[,i])
}


str(var.df)
tb<-var.df %>% group_by(groups) %>% summarize_all(funs(mean))
View(tb)

# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
bootvar.df<-t(var.df)
fit <- pvclust(bootvar.df, method.hclust="ward",method.dist="euclidean", nboot=1000,parallel=TRUE)
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)

seplot(fit)
seplot(fit, identify=TRUE)
print(result, which=c(21, 31))


# #################################################################
# # merge model evaluation csvs
# #################################################################

setwd(paste0(dir_out,"/ensemble-evals/"))
emval.files <- list.files(pattern="*EM_eval*")
emval.paths<-paste0(getwd(),"/",emval.files)

# ### merge variable importance (ensemble) csvs into dataframe
merged.em.evals <-
  do.call(rbind,
          lapply(emval.paths, read.csv))

View(merged.em.evals)
colnames(merged.em.evals)


### group variety and variable
dply<-merged.em.evals %>% group_by(variety,Model.name,Eval.metric) %>% summarize(Testing.data=mean(Testing.data),Sensitivity=mean(Sensitivity),Specificity=mean(Specificity))
View(dply)

dply<-as.data.frame(dply)
dply$Model.name<-sub("_MAXENT.Phillips_mergedRun_mergedData","",dply$Model.name)

var.avg<-dply %>% group_by(variety,Eval.metric,Model.name) %>% summarize(Testing.data=mean(Testing.data),Sensitivity=mean(Sensitivity),Specificity=mean(Specificity))
View(var.avg)

# model scores tend to be similar for each wmean ensemble
# average by evaluation metric
var.avg2<-dply %>% group_by(variety,Eval.metric) %>% summarize(Testing.data=mean(Testing.data),Sensitivity=mean(Sensitivity),Specificity=mean(Specificity))
View(var.avg2)

write.csv(var.avg2,file=paste0(dir_out,"/ensemble-evals/avg_score_by_metric_by_variety_with_ss.csv"))

# get mean sensitivity and specificity by variety
meanss<-var.avg2 %>% group_by(variety) %>% summarize(mean(Sensitivity),mean(Specificity))
meanss<-as.data.frame(meanss)
mean(meanss$`mean(Sensitivity)`)
#  78.5941
mean(meanss$`mean(Specificity)`)
#  91.15314
write.csv(meanss,file=paste0(dir_out,"/ensemble-evals/mean-spec-sens-by-variety.csv"))

spread(var.avg,Eval.metric,Testing.data)

### spread to get columns for each variety
colnames(var.avg)
var.avg2<-var.avg[,1:3]
tidy<-var.avg %>% group_by(variety,Eval.metric) %>% summarize(mean(Sensitivity),mean(Specificity)) %>% spread(key=c(Eval.metric),value=Testing.data)
View(tidy)
auc<-read.csv("E:\\thesis\\03_output\\ensemble-evals\\mean_auc_by_variety.csv")
tidy<-as.data.frame(tidy)
tidyauc<-merge(tidy,auc,by="variety")
mean(tidyauc$Mean_AUC)
# mean AUC 0.9392401
write.csv(tidyauc,file=paste0(dir_out,"/ensemble-evals/avg_score_by_variety_spread_by_metric-w-auc.csv"))

# #################################################################
# # merge range change dataframe
# #################################################################

setwd(paste0(dir_out,"/range-size/dfs"))
chg.files <- list.files(pattern="*")
chg.paths<-paste0(getwd(),"/",chg.files)

# ### merge variable importance (ensemble) csvs into dataframe
merged.rangesize <-
  do.call(rbind,
          lapply(chg.paths, read.csv))
library(tidyr)
library(dplyr)
colnames(merged.rangesize)
View(merged.rangesize)
colnames(merged.rangesize)

rangesorted <- merged.rangesize[order(merged.rangesize$time, merged.rangesize$species),]### group variety and variable


View(rangesorted)
colnames(rangesorted)


percchang<-rangesorted[,c("species","time","SpeciesRangeChange")]
percchang %>% select(SpeciesRangeChange,species,time)%>% spread(key=c(time),value=SpeciesRangeChange, convert=TRUE)
View(percchang)

colnames(percchang)[1]<-"variety"
colnames(percchang)[2]<-"Range_Change_%_(2041-2060)"
colnames(percchang)[3]<-"Range_Change_%_(2061-2070)"

write.csv(percchang,file=paste0(dir_range,"/range-change-by-period.csv"))
View(percchang2050)
# 
# 
# #################################################################
# # diversity zonal stats
# #################################################################
# #import required libraries
# library(maptools)
# library(raster)
# library(rgdal)
# par(mar = rep(2, 4))
# 
# dir_range<-paste0(dir_out,"/range-size")
# #list files (in this case raster TIFFs)
# alpha.div <- raster(paste0(dir_range,"/proj_current_alpha_diversity.grd"))
# alpha.div50 <- raster(paste0(dir_range,"/proj_rcp85_50_alpha_diversity.grd"))
# alpha.div70 <- raster(paste0(dir_range,"/proj_rcp85_70_alpha_diversity.grd"))
# 
# 
# #read-in the polygon shapefile
# poly <- readOGR("E:\\thesis\\01_data\\ind\\presindigw\\presindigw.shp",layer="presindigw")
# ob <- SpatialPolygons(poly@polygons,proj4string=poly@proj4string)
# ob<-spTransform(ob,proj4string(alpha.div))
# poly<-spTransform(poly,proj4string(alpha.div))
# 
# #extract raster cell count (sum) within each polygon area (poly)
# library(spatialEco)
# 
# mean.stat <- function(x, na.rm=TRUE) { 
#    if (na.rm) 
#        x <- x[!is.na(x)]
#   sum(x)/length(x) 
#   }  
# 
# med.stat <- function(x, na.rm=TRUE) { 
#    if (na.rm) 
#        x <- x[!is.na(x)]
# 	median(x)
#   }  
# 
# min.stat <- function(x, na.rm=TRUE) { 
#    if (na.rm) 
#        x <- x[!is.na(x)]
# 	min(x)
#   }  
# 
# df<-poly@data
# df<-as.data.frame(df)
# mean.div<-zonal.stats(poly, alpha.div, stat=mean.stat, trace = TRUE, plot = FALSE)
# med.div<-zonal.stats(poly, alpha.div, stat=med.stat, trace = TRUE, plot = FALSE)
# min.div<-zonal.stats(poly, alpha.div, stat=min.stat, trace = TRUE, plot = FALSE)
# 
# mean.div50<-zonal.stats(poly, alpha.div50, stat=mean.stat, trace = TRUE, plot = FALSE)
# med.div50<-zonal.stats(poly, alpha.div50, stat=med.stat, trace = TRUE, plot = FALSE)
# min.div50<-zonal.stats(poly, alpha.div50, stat=min.stat, trace = TRUE, plot = FALSE)
# 
# mean.div70<-zonal.stats(poly, alpha.div70, stat=mean.stat, trace = TRUE, plot = FALSE)
# med.div70<-zonal.stats(poly, alpha.div70, stat=med.stat, trace = TRUE, plot = FALSE)
# min.div70<-zonal.stats(poly, alpha.div70, stat=min.stat, trace = TRUE, plot = FALSE)
# 
# 
# length(mean.div)
# df$mean.div<-mean.div
# df$med.div<-med.div
# df$min.div<-min.div
# 
# df$mean.div50<-mean.div50
# df$med.div50<-med.div50
# df$min.div50<-min.div50
# 
# df$mean.div70<-mean.div70
# df$med.div70<-med.div70
# df$min.div70<-min.div70
# 
# poly<-readOGR("E:\\thesis\\01_data", layer="municipios-div-ind")
# df<-poly@data
# 
# 
# 
# colnames(df)
# df$mean.div50chg<-(df$men_dv50-df$mean_dv)/df$mean_dv
# 
# df<-df[complete.cases(df), ]
# 
# df$Presenc<-iconv(df$Presenc, from="UTF-8", to="LATIN1")
# 
# class(df$Presenc)
# df$Presenc<-factor(df$Presenc)
# df$Presenc
# df$Presenc<-factor(df$Presenc,levels=c("Sin poblaci󮠩ndna","Poblaci󮠩ndna dispersa","Poblaci󮠣on presencia indna","Poblaci󮠩ndna"),order=TRUE)
# 
# 
# 
# model = lm(mean.div50chg~ df$Presenc, 
#            data=df)
# 
# model 
# 
# library(car)
# 
# aov<-aov(model, Type="II",
#       white.adjust=TRUE)
# aov
# Anova<-Anova(model, Type="II",
#       white.adjust=TRUE)
# Anova
# 
# 
# thsd<-TukeyHSD(aov,ordered = TRUE, conf.level = 0.95)
# thsd
# # get categories of indigenous population levels
# df$indpobpct<- df$pobindi/df$pobtot
# df$indpobpct<-as.numeric(df$indpobpct)
# 
# df$ind_pct_group <- cut(df$indpobpct, 5)
# 
# str(df)
# nrow(poly@data)
# nrow(df)
# 
# poly2<-poly
# poly2@data<-merge(poly@data,df)
# nrow(poly2@data)
# length(unique(poly2@data$MUN_OFICIA))
# poly2@data$ID<-c(1:nrow(poly2@data))
# colnames(poly2@data)
# 
# 
# colnames(df)
# colnames(df)
# as.numeric(df$mean.div)
# df<-as.data.frame(df)
# View(df)
# nrow(df)
# 
# save.image(file="E:\\thesis\\01_data\\div_aov.RData")
# 
# 
# 
# 
# 
# poinmun <- readOGR("E:\\thesis\\01_data\\ind\\poinmun10gw\\poinmun10gw.shp",layer="poinmun10gw")
# 
# df.2<-merge(poinmun@data,df,by="MUN_OFICIA")
# 
# 
# View(df.2)
# 
# #write to a CSV file
# names(df.2) <- tolower(names(df.2))
# 
# write.csv(df.2, file = paste0(dir_out,"/alpha-div-reg-df.csv"))
# str(df.2)
# names(df.2) <- tolower(names(df.2))
# library(leaps)
# 
# 
# regsubsets.out <-
#     regsubsets(min.div ~ p3t10+ p3i10+ mon10+ bil10+ nei10+  nli10+ neli10+area.x+perimeter.x+presencia + pobindi+marmun ,
#                data = df.2,
#                nbest = 1,       # 1 best model for each number of predictors
#                nvmax = NULL,    # NULL for no limit on number of variables
#                force.in = NULL, force.out = NULL,
#                method = "exhaustive")
# regsubsets.out
# 
# summary.out <- summary(regsubsets.out)
# as.data.frame(summary.out$outmat)
# 
# plot(regsubsets.out, scale = "adjr2", main = "Adjusted R^2")
# library(car)
# layout(matrix(1:2, ncol = 2))
# ## Adjusted R2
# res.legend <-
#     subsets(regsubsets.out, statistic="adjr2", legend = FALSE, min.size = 5, main = "Adjusted R^2")
# ## Mallow Cp
# res.legend <-
#     subsets(regsubsets.out, statistic="cp", legend = FALSE, min.size = 5, main = "Mallow Cp")
# 
# which.max(summary.out$adjr2)
# 
# best.model <- lm(mea.div ~ p3t10+ p3i10+  nei10+  neli10+area.x+perimeter.x+presencia*marmun+ pobindi , data = df.2)
# summary(best.model)
# 
# 
# summary.out$which[11,]
# abline(a = 1, b = 1, lty = 2)
# library(corrgram)
# library("ggpubr")
# ggscatter(df, x = "index_1", y = "pobpct", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "pearson",
#           xlab = "Maize Alpha-diversity", ylab = "Indigenous Population (%)")
# 
# 
# # Shapiro-Wilk normality tests 
# shapiro.test(df$index_1) # => p = 0.1229
# 
# shapiro.test(df$pobpct) # => p = 0.09
# 
# library("ggpubr")
# # mpg
# ggqqplot(df$index_1, ylab = "Alpha Diversity")
# # wt
# ggqqplot(df$pobpct, ylab = "Ind. Pob. (%)")
# 
# cor.pearson <-cor.test(df$index_1,df$pobpct,  method = c("pearson"))
# cor.pearson
# cor.kendall <-cor.test(df$index_1,df$pobpct,  method = "kendall")
# cor.kendall
# cor.spearman <-cor.test(df$index_1,df$pobpct,  method = "spearman")
# cor.spearman
# 
# 
# 
# # Figner-Killeen Test of Homogeneity of Variances
# colnames(df)
# str(df)
# 
# head(df$group)
# fligner.test(index_1~group, data=df)
# bartlett.test(index_1~group, data = df)
# kruskal.test(index_1~group, data = df)
# View(df$index_1)
# 
# 
# 
# plot(alpha.div)
# pobindras<-raster("E:\\thesis\\01_data\\ind\\pob-ind.grd")
# library(spatialEco)
# pobindras<-projectRaster(pobindras,alpha.div)
# pobindras<-resample(pobindras,alpha.div,method="bilinear")
# 
# rascorr3<-rasterCorrelation(pobindras, alpha.div, s = 3, type = "pearson", file.name = NULL)
# rascorr3<-rascorr
# plot(rascorr3)
# rascorr5<-rasterCorrelation(pobindras, alpha.div, s = 5, type = "pearson", file.name = NULL)
# rascorr7<-rasterCorrelation(pobindras, alpha.div, s = 7, type = "pearson", file.name = NULL)
# rascorr25<-rasterCorrelation(pobindras, alpha.div, s = 25, type = "pearson", file.name = NULL)
# 
# 
# library(raster)
# library(SDMTools)
# library(adehabitat)
# 
# rAsc <- asc.from.raster(pobindras) # Function from SDMTools to convert to asc format
# pobindrasspdf <- asc2spixdf(rAsc)
# 
# rAsc <- asc.from.raster(alpha.div) # Function from SDMTools to convert to asc format
# alpha.divspdf<- asc2spixdf(rAsc)
# 
# class(alpha.div)
# modttest<-spatialEco::raster.modifed.ttest(pobindrasspdf , alpha.divspdf, d = "AUTO",
#   sub.sample = TRUE, type = "hexagon", p = 0.1)
# 
# modttest
