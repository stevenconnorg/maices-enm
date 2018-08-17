source("init.R")

#################################################################
# PRELIMINARY DATA FORMATTING
#################################################################

# read in model raster stacks in LONG/LAT
# can't project onto equal area because Maxent needs resolution with equal x, y
# this data came from lat/long, so the projected 
library(raster)
library(rgdal)
library(biomod2)


setwd(dir_bm) 
#setwd(dir_bmf) 


# get species vector name by getting all with ensemble projections for 2070 time period
# or whichever was last in the BioModApply function/loop

emprojfiles<-Sys.glob(paste0(getwd(),"/*/*/","*70.ensemble.projection.out")) # get ensemble files out
emproj.dirs<-sub(paste0(getwd(),"/"),"",emprojfiles) # remove everything from working directory path
sp.n<-sub(" */.*", "",emproj.dirs) # remove everything before first slash to get variety names that have ensemble models
sp.n


#################################################################
# INITIALIZE FUNCTION TO APPLY TO EACH VARIETY
#################################################################

# dir_out has data from the target-specific stuff
# dir_out<-dir_outf


dir_wmean<-paste0(dir_out,"/wmean-consensus")
dir_binary<-paste0(dir_out,"/binaries")
dir_rangeras<-paste0(dir_out,"/range-size/rasters")
dir_range<-paste0(dir_out,"/range-size")
dir_rangedf<-paste0(dir_out,"/range-size/dfs")
dir_rangeR<-paste0(dir_out,"/range-size/RData")
dir_div<-paste0(dir_out,"/diversity")

fun.auc.df<-data.frame(Landrace = as.character(),
                      AUC =as.numeric(),
                       x.Sens =as.numeric(),
                       x.Spec =as.numeric(),
                       SS_min_dif =as.numeric(),
                       SS_max_sum =as.numeric(),
                       Min_Err =as.numeric(),
                       Min_Err_Cut=as.numeric())

#PostBioEval <-function(sp) {
for (sp in sp.n){
  #tryCatch({

  setwd(dir_bm)
  modelfile<-load(paste0(getwd(),"/",sp,"/",sp,".",sp,"_current.models.out"))
  myBiomodModelOut<-get(modelfile)
  

  emmodelfile<-load(paste0(getwd(),"/",sp,"/",sp,".",sp,"_currentensemble.models.out"))
  myBiomodEM<-get(emmodelfile)
  
  proj_current<-load(paste0(getwd(),"/",sp,"/proj_current/",sp,".current.ensemble.projection.out"))
  myBiomodProj<-get(paste0(sp,".current.ensemble.projection.out"))
  myBiomod_raster <- get_predictions(myBiomodProj)
  
  
  KAPPAem<-paste0(sp,"_EMwmeanByKAPPA_MAXENT.Phillips_mergedRun_mergedData")
  TSSem<-paste0(sp,"_EMwmeanByTSS_MAXENT.Phillips_mergedRun_mergedData")
  ROCem<-paste0(sp,"_EMwmeanByROC_MAXENT.Phillips_mergedRun_mergedData")
  
  
  # # # if working on local laptop, you have to redirect file name
  rasKAPPA<-subset(myBiomod_raster,KAPPAem)
  slotname<-slot(rasKAPPA@file,"name")
  newname<-sub("/gpfs/home/scg67/thesis/",paste0(root,"/"),slotname)
  slot(rasKAPPA@file,"name")<-newname
  rasKAPPA
  
  rasTSS<-subset(myBiomod_raster,TSSem)
  slotname<-slot(rasTSS@file,"name")
  newname<-sub("/gpfs/home/scg67/thesis/",paste0(root,"/"),slotname)
  slot(rasTSS@file,"name")<-newname
  rasTSS
  
  rasROC<-subset(myBiomod_raster,ROCem)
  slotname<-slot(rasROC@file,"name")
  newname<-sub("/gpfs/home/scg67/thesis/",paste0(root,"/"),slotname)
  slot(rasROC@file,"name")<-newname
  rasROC
  
  
  pa<-data.frame(read.csv(file=paste0(dir_out,"/pa_dataframe_EA.csv")))
  colnames(pa)
  colnames(pa) <- sub("_", ".", colnames(pa))
  colnames(pa) <- sub("_", ".", colnames(pa))
  colnames(pa) <- sub("_", ".", colnames(pa))
  colnames(pa)
  
  
  myBiomod_raster<-stack(rasKAPPA,rasTSS,rasROC)
  points<-SpatialPoints(coords = pa[,2:3])
  ex<-raster::extract(myBiomod_raster,points)
  ex<-ex/1000
  colnames(ex)

  myResp <- as.numeric(pa[,sp])
  spp_occ <- myResp
  library(ROCR)
  ROC_curve <- data.frame(result=ex, spp=spp_occ) #préparation
  ROC_curve[is.na(ROC_curve)] <- 0
  colnames(ROC_curve)
  ROC_curve<-ROC_curve[complete.cases(ROC_curve),]
  colnames(ROC_curve)
  ROC_curve[, !colnames(ROC_curve) %in% KAPPA]
  KAPPA<-paste0("result.",sp,"_EMwmeanByKAPPA_MAXENT.Phillips_mergedRun_mergedData")
  TSS<-paste0("result.",sp,"_EMwmeanByTSS_MAXENT.Phillips_mergedRun_mergedData")
  ROC<-paste0("result.",sp,"_EMwmeanByROC_MAXENT.Phillips_mergedRun_mergedData")
  nrow(ROC_curve[, !colnames(ROC_curve) %in% KAPPA])
  length(ROC_curve$spp)
  class(unlist(ROC_curve[, !colnames(ROC_curve) %in% KAPPA]))
  class(ROC_curve[, !colnames(ROC_curve) %in% KAPPA])
  
  
  
  # ROC_curve_KAPPA<- ROCR::prediction(ROC_curve[, !colnames(ROC_curve) %in% KAPPA][,2], as.vector(ROC_curve$spp)) #run
  # ROC_curve_TSS <- prediction(ROC_curve[, !colnames(ROC_curve) %in% TSS][,2], ROC_curve$spp) #run
  ROC_curve_ROC <- prediction(ROC_curve[, !colnames(ROC_curve) %in% ROC][,2], ROC_curve$spp) #run
  
   dir_rocr<-paste0(dir_out,"/ensemble-evals/rocr/RData")
  dir.create(dir_rocr,recursive=TRUE,showWarnings=FALSE)
  save(ROC_curve_ROC,file=paste0(dir_rocr,"/",sp,"-ROC-ens-ROC-curve-prediction.RData"))
  
  
  # https://davidrroberts.wordpress.com/2015/09/22/quick-auc-function-in-r-with-rocr-package/
  # AUC function
  fun.auc <- function(pred,obs){
    # Run the ROCR functions for AUC calculation
    ROC_perf <- performance(prediction(pred,obs),"tpr","fpr")
    ROC_sens <- performance(prediction(pred,obs),"sens","spec")
    ROC_err <- performance(prediction(pred, labels=obs),"err")
    ROC_auc <- performance(prediction(pred,obs),"auc")
    # AUC value
    AUC <- ROC_auc@y.values[[1]] # AUC
    # Mean sensitivity across all cutoffs
    x.Sens <- mean(as.data.frame(ROC_sens@y.values)[,1])
    # Mean specificity across all cutoffs
    x.Spec <- mean(as.data.frame(ROC_sens@x.values)[,1])
    # Sens-Spec table to estimate threshold cutoffs
    SS <- data.frame(SENS=as.data.frame(ROC_sens@y.values)[,1],SPEC=as.data.frame(ROC_sens@x.values)[,1])
    # Threshold cutoff with min difference between Sens and Spec
    SS_min_dif <- ROC_perf@alpha.values[[1]][which.min(abs(SS$SENS-SS$SPEC))]
    # Threshold cutoff with max sum of Sens and Spec
    SS_max_sum <- ROC_perf@alpha.values[[1]][which.max(rowSums(SS[c("SENS","SPEC")]))]
    # Min error rate
    Min_Err <- min(ROC_err@y.values[[1]])
    # Threshold cutoff resulting in min error rate
    Min_Err_Cut <- ROC_err@x.values[[1]][which(ROC_err@y.values[[1]]==Min_Err)][1]
    # Kick out the values
    round(cbind(AUC,x.Sens,x.Spec,SS_min_dif,SS_max_sum,Min_Err,Min_Err_Cut),3)

    
      # png(paste0(dir_rocr,"/AUC/",sp,"-AUC.png"), width =10, height = 13, units = 'in', res = 600)
      #
      # plot(ROC_perf,colorize=T,print.cutoffs.at=seq(0,1,by=0.1),lwd=3,las=1)
      # # Add some statistics to the plot
      # text(1,0.30,labels=paste(sp),adj=1)
      # text(1,0.25,labels=paste("Npres = ",sum(obs==1),sep=""),adj=1)
      # text(1,0.20,labels=paste("Nabs = ",sum(obs==0),sep=""),adj=1)
      # text(1,0.15,labels=paste("AUC = ",round(ROC_auc@y.values[[1]],digits=2),sep=""),adj=1)
      # text(1,0.10,labels=paste("Sens = ",round(mean(as.data.frame(ROC_sens@y.values)[,1]),digits=2),sep=""),adj=1)
      # text(1,0.05,labels=paste("Spec = ",round(mean(as.data.frame(ROC_sens@x.values)[,1]),digits=2),sep=""),adj=1)
      # abline(a=0, b= 1)
      # dev.off()
      
      
      
      
      
      # dir.create(paste0(dir_rocr,"/AUC/"),recursive=TRUE,showWarnings=FALSE)
      # dir.create(paste0(dir_rocr,"/Accuracy/"),recursive=TRUE,showWarnings=FALSE)
      #
    
  }
  
  

  
  
 #  Run the function with the example data
   
  pred<-ROC_curve[, !colnames(ROC_curve) %in% ROC][,2]
  obs<- as.vector(ROC_curve$spp)
  # # #
  # # 
  # # # bind stats to dataframe
  auc.dat<-fun.auc(pred, obs)
  fun.auc.df<-rbind(fun.auc.df,data.frame(sp,auc.dat))
  # #
  




}

write.csv(fun.auc.df,paste0(dir_rocr,"/AUC_stats.csv"))



sp<-sp.n[1]
for (sp in sp.n){
  dir_rocrDat<-paste0(dir_out,"/ensemble-evals/rocr/RData")
  
  load(paste0(dir_rocrDat,"/",sp,"-ROC-ens-ROC-curve-prediction.RData"))
  
  # get AUC plots
  dir_rocr<-paste0(dir_out,"/ensemble-evals/rocr")
  
  dir.create(paste0(dir_rocr,"/Accuracy/"),recursive=TRUE,showWarnings=FALSE)
   perf <- performance(ROC_curve_ROC, "tpr", "fpr")

  png(paste0(dir_rocr,"/Accuracy/",sp,"-Accuracy.png"), width =10, height = 13, units = 'in', res = 600)

   plot(perf)
   abline(a=0, b= 1)
   dev.off()
}

  # KAPPA_acc <- performance( ROC_curve_KAPPA, "acc" )@y.values
  # TSS_acc <- performance( ROC_curve_TSS, "acc" )@y.values
  # ROC_acc <- performance( ROC_curve_ROC, "acc" )@y.values
  # 
  # KAPPA_err<- performance( ROC_curve_KAPPA, "err" )@y.values
  # TSS_err <- performance( ROC_curve_TSS, "err" )@y.values
  # ROC_err <- performance( ROC_curve_ROC, "err" )@y.values
  # 
  # KAPPA_fpr<- performance( ROC_curve_KAPPA, "fpr" )@y.values
  # TSS_fpr <- performance( ROC_curve_TSS, "fpr" )@y.values
  # ROC_fpr <- performance( ROC_curve_ROC, "fpr" )@y.values
  # 
  # KAPPA_tpr<- performance( ROC_curve_KAPPA, "tpr" )@y.values
  # TSS_tpr <- performance( ROC_curve_TSS, "tpr" )@y.values
  # ROC_tpr <- performance( ROC_curve_ROC, "tpr" )@y.values
  # 
  # KAPPA_fnr<- performance( ROC_curve_KAPPA, "fnr" )@y.values
  # TSS_fnr <- performance( ROC_curve_TSS, "fnr" )@y.values
  # ROC_fnr <- performance( ROC_curve_ROC, "fnr" )@y.values
  # 
  # KAPPA_f<- performance( ROC_curve_KAPPA, "f" )@y.values
  # TSS_f <- performance( ROC_curve_TSS, "f" )@y.values
  # ROC_f <- performance( ROC_curve_ROC, "f" )@y.values
  # 
  # KAPPA_sar<- performance( ROC_curve_KAPPA, "sar" )@y.values
  # TSS_sar <- performance( ROC_curve_TSS, "sar" )@y.values
  # ROC_sar <- performance( ROC_curve_ROC, "sar" )@y.values
  # 

  eval.df<-data.frame(Model=unique(enevalmods$Model.name),
                      AUC=rbind(KAPPA_AUC,TSS_AUC,ROC_AUC
                      )
                      
  )
  eval.df$variety<-sp
  View(eval.df)
  rownames(AUC.df)<-1:nrow(AUC.df)
  class(AUC.df)
  AUC.df[,2]<-as.numeric(AUC.df[,2])
  write.csv(AUC.df,file=paste0(dir_out,"/ensemble-evals/",sp,"-ensemble-AUC.csv"))
  