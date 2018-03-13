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


PostBioEval <-function(sp) {
  #tryCatch({

  #################################################################
  # CAPTURE MODEL RESPONSE CURVES TO FILE BY ALGORITHM
  #################################################################
  

  #
  setwd(dir_bm)
  modelfile<-load(paste0(getwd(),"/",sp,"/",sp,".",sp,"_current.models.out"))
  myBiomodModelOut<-get(modelfile)

  #################################################################
  # CAPTURE DATA INPUT AND MODEL OUTPUTS
  ################################################################
  
  # write data used for modelling
  capture.output(get_formal_data(myBiomodModelOut),
                 file=paste0(dir_out,"/",myRespName,"/",myRespName,"_model_data.txt"))
  
  
  #################################################################
  # GET MODEL EVALUATIONS
  #################################################################
  
  #print(paste0("Capturing Model Scores Graph for ",sp))
  dir_evals<-paste0(dir_out,"/model-evals/")
  dir.create(dir_evals,recursive=TRUE,showWarnings=FALSE)
  
  
  # write em models built
  capture.output(get_built_models (myBiomodModelOut),
                 file=paste0(dir_evals,"/",sp,"_models.txt"))
  
  # capture em model evals
  print(paste0("Capturing Models Evaluations ",sp))
  
  
  evalmods<-get_evaluations(myBiomodModelOut,as.data.frame=TRUE)
  evalmods$variety<-sp
  write.csv(evalmods,file=paste0(dir_evals,"/",sp,"_models_eval.csv"))
  

  ### eval current model
  
  #print(paste0("Capturing Model Ensemble Evaluations for ",sp))
  evalmods<-get_evaluations(myBiomodModelOut,as.data.frame=TRUE)
  write.csv(evalmods,file=paste0(dir_evals,"/",sp,"_evals-df.csv"))
  
  #################################################################
  # GET ENSEMBLE MODEL EVALUATIONS
  #################################################################
  # 
  emmodelfile<-load(paste0(getwd(),"/",sp,"/",sp,".",sp,"_currentensemble.models.out"))
  myBiomodEM<-get(emmodelfile)
  #print(paste0("Capturing Model Scores Graph for ",sp))
  dir_enevals<-paste0(dir_out,"/ensemble-evals/")
  dir.create(dir_enevals,recursive=TRUE,showWarnings=FALSE)
  
  
  print(paste0("Capturing Ensemble Model Outputs for  ",sp))
  
  # write em models built
  capture.output(get_built_models (myBiomodEM),
                 file=paste0(dir_enevals,"/",sp,"_em_models.txt"))
  
  # capture em model evals
  
  print(paste0("Capturing Ensemble Models Evaluations ",sp))
  evalmodsEM<-get_evaluations(myBiomodEM,as.data.frame=TRUE)
  evalmodsEM$variety<-sp
  write.csv(evalmodsEM,file=paste0(dir_enevals,"/",sp,"_models_EM_eval.csv"))
  

  #################################################################
  # CAPTURE MESS STACKS
  #################################################################
  setwd(dir_bm)
  modelfile<-load(paste0(getwd(),"/",sp,"/",sp,".",sp,"_current.models.out"))
  myBiomodModelOut<-get(modelfile)
  
  # dat<-load(slot(myBiomodModelOut@formated.input.data,"link"))
  # coords<-data@coord
  # coordinates<-coordinates(coords)
  # 
  # 
  pa<-data.frame(read.csv(file=paste0(dir_out,"/pa_dataframe_EA.csv")))
  colnames(pa)
  colnames(pa) <- sub("_", ".", colnames(pa))
  colnames(pa) <- sub("_", ".", colnames(pa))
  colnames(pa) <- sub("_", ".", colnames(pa))
  colnames(pa)
  
  resp.occ.id <- which(pa[, sp] == 1)
  reference_points <- pa[resp.occ.id, c("Longitude.x.","Latitude.y.")]
  coords<-coordinates(reference_points)
  

# for each variety  
    grdnames<-list("present_modstack_EA1km","f50_modstack_EA1km","f70_modstack_EA1km")
    
    for (i in grdnames){
      stack<-stack(paste0(dir_stacks,i,".grd"))
      reference_points<-extract(stack,coords)
      dir_MESS<-paste0(dir_out,"/",i,"-MESS")
      dir.create(dir_MESS,showWarnings=FALSE,recursive = T)
      
      mess.out <- mess(x=stack, v=reference_points, full=FALSE,paste0(dir_MESS,"/",sp,"-MESS.grd"))
    }

#       
    reference_points<-as.data.frame(stack)

    for (i in grdnames){
      stack<-stack(paste0(dir_stacks,i,".grd"))
      dir_MESS<-paste0(dir_out,"/",i,"-MESS")
      dir.create(dir_MESS,showWarnings=FALSE,recursive = T)
      
      mess.out <- mess(x=stack, v=reference_points, full=FALSE,paste0(dir_MESS,"/",sp,"-MESS.grd"))
    }




    grdnames<-list("present_modstack_EA1km","f50_modstack_EA1km","f70_modstack_EA1km")
    
    pstack<-stack(paste0(dir_stacks,grdnames[1],".grd"))[[1:10]]
    names(pstack)
    f50stack<-stack(paste0(dir_stacks,grdnames[2],".grd"))[[1:10]]
    
    f70stack<-stack(paste0(dir_stacks,grdnames[3],".grd"))[[1:10]]
    
      refstack<-stack(paste0(dir_stacks,grdnames[1],".grd"))[[1:10]]
      
      xy<-xyFromCell(refstack)
      
      
      
      # for each project (full target + background dataset)
      
      latlong <- pa[, c("Longitude.x.","Latitude.y.")]
      coords<-coordinates(latlong)
      pstack<-stack(paste0(dir_stacks,grdnames[1],".grd"))[[1:10]]
      
      reference_points<-extract(coords,pstack)
      
    for (i in grdnames){
      reference_points<- refstack
      
      stack<-stack(paste0(dir_stacks,i,".grd"))[[1:10]]

      dir_MESS<-paste0(dir_out,"/MESS/",i,"-MESS")
      dir.create(dir_MESS,showWarnings=FALSE,recursive = T)
      mess.out <- mess(x=stack, v=reference_points, full=TRUE,paste0(dir_MESS,"/bioclimatic-full-raster-MESS.grd"))
	  
  }
  
  
  # get average mess index across each period

MESS_stats<-data.frame(Present_MESS=as.numeric(),f2050_MESS=as.numeric(),f2070_MESS=as.numeric())
    grdnames<-list("present_modstack_EA1km","f50_modstack_EA1km","f70_modstack_EA1km")
    for (i in grdnames){
      dir_MESS<-paste0(dir_out,"/",i,"-MESS")
      files<-list.files(dir_MESS,pattern="full-raster-MESS.grd",full.names=TRUE)
	  mess_stack<-stack(files)
	  names(stack)<-sp.n
	  mean<-as.data.frame(cellStats(mess_stack,'mean'))
	  colnames(mean)[1]<-paste0('mean_MESS',i)
	  rownames(mean)<-sp.n
	  write.csv(mean,paste0(dir_MESS,"/average-MESS-across-landraces.csv"))
	  #meanMess<-stackApply(mess_stack,indices=rep(1,nlayers(mess_stack)),fun=mean,na.rm=TRUE)
	  #cellStats(meanMess,stat='mean')
  }

  
    
    MESSfiles<-Sys.glob(paste0(root,"/03_output/MESS/*/*target*.grd"))
    
    
    files<-list.files(paste0(root,"/03_output/MESS/"),pattern="full-raster-MESS.grd",full.names=TRUE,recursive=TRUE)
    files<-list.files(paste0(root,"/03_output/MESS/"),pattern="bioclimatic-full-raster-MESS.grd",full.names=TRUE,recursive=TRUE)
    
    messStack<-stack()
    for (mess in MESSfiles){
      messfile<-mess
      mess<-stack(messfile)
      mess<-mess
      levelplot(mess)
      messplt<-levelplot(s2[[c(3,2)]], names.attr=c("1970-2000","2061-2080"), margin = list(FUN = 'median'),
                         xlab="X", ylab="Y",main=paste0("Biocliimatic Multivariate Environmental Similarity Surfaces"))
      
      messplt<-levelplot(mess, names.attr=sub("_mess","",names(mess)), margin = list(FUN = 'median'),
                         xlab="X", ylab="Y",main=paste0("Multivariate Environmental Similarity Surfaces"))
      
      png(paste0(dir_out,"/MESS/MESS-bio-levelplot.png"), width = 10, height = 10, units = 'in', res = 600)
      print(messplt)
      dev.off()
      
      # messStack<-stack(messStack,mess)
    }
    library(rasterVis)
    names(messStack)

    messplt<-levelplot(messStack, names.attr=c("1970-2000","2041-2060","2061-2080"), margin = list(FUN = 'median'),
              xlab="X", ylab="Y",main=paste0("Multivariate Environmental Similarity Surfaces"))
    
    png(paste0(dir_out,"/MESS/MESS-levelplot.png"), width = 10, height = 10, units = 'in', res = 600)
    print(messplt)
    dev.off()
    
  #################################################################
  # GET VARIABLE IMPORTANCE
  #################################################################
  # 
  # 
  # print(paste0("Capturing Model Scores Graph for ",sp))
   dir_var<-paste0(dir_out,"/var-importance/")

   dir.create(dir_var,recursive=TRUE,showWarnings=FALSE)

   
  modvarimp<-get_variables_importance(myBiomodModelOut)
  modvarimp.df<-as.data.frame(modvarimp)
  modvarimp.df$variety<-sp
  write.csv(modvarimp.df,file=paste0(dir_var,"/",sp,"_model-var-imp-df.csv"))

  EMvarimp<-get_variables_importance(myBiomodEM)
  EMvarimp.df<-as.data.frame(rowMeans(EMvarimp))
  EMvarimp.df$variety<-sp

  write.csv(EMvarimp.df,file=paste0(dir_var,"/",sp,"_EM-mean-var-imp-df.csv"))

  #################################################################
  # CAPTURE MODEL RESPONSE CURVES TO FILE BY MODEL ALGORITHM RUN
  #################################################################

    print(paste0("Capturing Model Scores Graph for ",sp))
    dir_curves<-paste0(dir_out,"/response-curves/")
    dir.create(dir_curves,recursive=TRUE,showWarnings=FALSE)

    modelsubset<-BIOMOD_LoadModels(myBiomodModelOut,models=mod)

    pdf(width = 6, height = 6)
    response.plot2(models=modelsubset,
                   Data = get_formal_data(myBiomodModelOut,'expl.var'),
                   show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
                   save.file = "pdf",
                   col=rainbow(length(myBiomodModelOut@models.computed)),
                   name=paste0(dir_curves,sp,"_",mod,"_curves"),
                   do.bivariate=FALSE,
                   fixed.var.metric = 'median',
                   data_species = get_formal_data(myBiomodModelOut,'resp.var'),
                   ImageSize=5000,
                   plot=TRUE)
  
  
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
  #name=paste0(dir_curves,sp,"_",mod,"_",i,"_3d_curves"),
  #ImageSize=500,
  # plot=TRUE)
  
  
  # } # response plot3d loop by individual model
  
  
  
  #################################################################
  # GET MODEL SCORES GRAPH   -- under construction
  #################################################################
  
  
  
  #print(paste0("Capturing Model Scores Graph for ",sp))
  #dir_modscores<-paste0(dir_out,"/model-score-graphs/")
  #dir.create(dir_modscores,recursive=TRUE,showWarnings=FALSE)
  
  
  # for (metrics in seq_along(metrics)) ....
  
  #pdf(filename=paste0(dir_modscores,"/",sp,"_model_scores-",paste(metrics[1:2], collapse = "-"),".pdf"))
  #models_scores_graph(myBiomodModelOut,metrics = metrics[1:2],plot = TRUE)
  #dev.off()
  
  
  #pdf(filename=paste0(dir_modscores,"/",sp,"_model_scores-",paste(metrics[2:3], collapse = "-"),".pdf"))
  #models_scores_graph(myBiomodModelOut,metrics = metrics[2:3],by = 'cv_run',plot = TRUE, width = 12, height = 17, family = "Helvetica")
  #dev.off()
  
  
  #pdf(filename=paste0(dir_modscores,"/",sp,"_model_scores-",paste(metrics[c(1,3)], collapse = "-"),".pdf"))
  #models_scores_graph(myBiomodModelOut,metrics = metrics[c(1,3)],by = 'cv_run',plot = TRUE, width = 12, height = 17, family = "Helvetica")
  #dev.off()
  
  library(rasterVis)
  library(sp)
}



PostBioApply <-function(sp) {
  for (pr in projects){
    dir_wmean<-paste0(dir_out,"/wmean-consensus")
    dir.create(dir_wmean,recursive=TRUE,showWarnings=FALSE)
    
    #################################################################
    # GET AVERAGE WEIGHTED MEAN PREDICTIONS
    #################################################################
    
    # average weighted mean ensemble projections

    ## if (pr == projects[2:3]){
    # dir_bm<-dir_bmf
    #return(dir_bm)
    # } 
    
      emstack<-raster::stack(paste0(dir_bm,"/",sp,"/",pr,"/",pr,"_",sp,"_ensemble.grd"))
      return(emstack)

    emwmeans <- raster::subset(emstack, grep('EMwmean', names(emstack), value = T))
    
    
    emwmean<-raster::stackApply(emwmeans,indices=c(rep(1,raster::nlayers(emwmeans))),fun=mean)
    raster::writeRaster(emwmean,file=paste0(dir_wmean,"/",sp,"_",pr,"_em-wmean.grd"),format="raster",overwrite=T)
    
    
    #################################################################
    # GET BINARY PREDICTIONS FROM AVERAGE WEIGHTED MEAN ENSEMBLES
    #################################################################
    
    # matrix to reclassify binary total consensus averages
    class.m <- c(0, 750, 0,
                 750,1001, 1,
                 NA, NA,  NA)
    
    # reshape the object into a matrix with columns and rows
    rcl.m <- matrix(class.m, ncol=3, byrow=TRUE)
    
    
    wmean_bin <- reclassify(emwmean, 
                            rcl.m)
    
    raster::writeRaster(wmean_bin,file=paste0(dir_binary,"/binary_750_wmean_avg_",sp,"_",pr,".grd"),overwrite=TRUE)
    
  }
  
  #################################################################
  # GET PLOT OF CURRENT WEIGHTED MEAN PREDICTION
  #################################################################
  library(rasterVis)
  wmean<-Sys.glob(paste0(dir_wmean,"/",sp,"*.grd")) # get ensemble files out
  current.wmean<-raster(wmean[1])
  mapTheme <- rasterTheme(region=brewer.pal(5,"Greens"))
  
  
  plt <-levelplot(current.wmean, margin = list(FUN = 'median'), par.settings=mapTheme,
                  xlab="X", ylab="Y",main=paste0(sp," : 1970 - 2000 \nProjected Relative Probability of Presence") )
  
  png(paste0(dir_wmean,"/",sp,"-current-wmean-em-avg.png"), width = 10, height = 10, units = 'in', res = 600)
  print(plt  )
  dev.off()
  
  f50.wmean<-raster(wmean[2])
  f70.wmean<-raster(wmean[3])
  
  f50.wmean<-resample(f50.wmean,current.wmean,method="bilinear")
  f70.wmean<-resample(f70.wmean,current.wmean,method="bilinear")
  
  wmstack<-stack(current.wmean,f50.wmean,f70.wmean)
  names(wmstack)<-c("1970-2000","2041-2060","2061-2080")
  wmstack<-wmstack/1000
 
  
  
  plt <-levelplot(wmstack, margin = list(FUN = 'median'), par.settings=mapTheme,
                  xlab="X", ylab="Y",main=paste0(sp," - \nProjected Relative Probability of Presence"), names.attr=c("1970-2000","2041-2060","2061-2080")
  )
  
  png(paste0(dir_wmean,"/",sp,"-wfuture-wmean-em-avg.png"), width = 10, height = 10, units = 'in', res = 600)
  print(plt)
  dev.off()
}

 # biomodpostapply
  

library(Rmpi)
library(parallel)
library(snowfall)
ncpus<-detectCores()
sfInit( parallel=TRUE, cpus=ncpus, type="MPI",
        slaveOutfile = paste0(Sys.getenv("SLURM_JOB_ID"),"_sfInit.log"))

# # Export packages to snowfall
print("Exporting packages to cluster")
sfLibrary('raster', character.only = TRUE)
sfLibrary('rgdal', character.only = TRUE)
sfLibrary('gtable', character.only = TRUE)
sfLibrary('grid', character.only = TRUE)
sfLibrary('ggplot2', character.only = TRUE)
sfLibrary('rasterVis', character.only = TRUE)

# # Export environment to snowfall
sfExportAll()

sfLapply(sp.n, PostBioApply)

sfStop()

#######

#################################################################
# GET TOTAL CONSENSUS ENSEMBLE BINARY PREDICTIONS
#################################################################

rasterOptions(maxmemory=3e+08, chunksize=1e+06)


for (p in projs){
    for (sp in sp.n){
      
        ## if (p == projs[2:3]){
        # dir_bm<-dir_bmf
        #return(dir_bm)
        # } 
      
         bin_stack<-stack(Sys.glob(paste(dir_bm,"/",sp,"/proj_",p,"/proj_",p,"_",sp,"_ensemble_*bin.grd",sep="")))
         #bin_stack<-brick(bin_stack)
         
         
         bin.tc<-raster::stackApply(bin_stack,indices=c(rep(1,raster::nlayers(bin_stack))),fun=mean,na.rm=TRUE)
                  # reclassify the raster using the reclass object - rcl.m
         
         
         
         # matrix to reclassify binary total consensus averages
         class.m <- c(0, .700, 0,
                      .701,1.1, 1,
                      NA, NA,  NA)
         
         # reshape the object into a matrix with columns and rows
         rcl.m <- matrix(class.m, ncol=3, byrow=TRUE)
         
         
         bin.consensus <- reclassify(bin.tc, 
                              rcl.m)
         writeRaster(bin.tc,file=paste0(dir_binary,"/binary_tc-avg_",sp,"_",p,".grd"),overwrite=TRUE)
         writeRaster(bin.consensus,file=paste0(dir_binary,"/binary_tc_70threshold",sp,"_",p,".grd"),overwrite=TRUE)
 
      } # sp in sp.n

} # for p in projects

        
#################################################################
# GET ALPHA DIVERSITY
# and get mean and sum alpha diversity by aggregation
#################################################################


for (p in projs){
  tmpl_ras<-raster(Sys.glob(paste0(dir_binary,"/binary_tc_70threshold*",p,"*.grd"))[1])
  
  divstack<-stack()
  divstack<-resample(divstack,tmpl_ras,method="bilinear")
  
  files<-Sys.glob(paste0(dir_binary,"/binary_tc_70threshold*",p,"*.grd"))
  
  for (f in files){
    ras<-raster(f)
    ras<-resample(ras,tmpl_ras,method="bilinear")
    divstack <- stack( divstack , ras )  
    }
  
  # divstack<-brick(Sys.glob(paste0(dir_binary,"/binary_tc_70threshold*",p,"*.grd")))
  # for (i in c(1:nlayers(divstack))){
  #   divstack[[i]]<-resample(divstack[[i]],tmpl_ras,method="bilinear")
  # }

  alphadiv<-raster::stackApply(divstack,indices=c(rep(1,raster::nlayers(divstack))),fun=sum)
  
    writeRaster(alphadiv,file=paste0(dir_div,"/maize_alpha_diversity_",p,".grd"),overwrite=TRUE)
  
  for (dim in c(2.5,5,7,10,15,20)){
    
      for (stat in c("mean","sum","min","max")){            
        area<-dim[1]*dim[1]
      area<-round(area,digits = 0)
      
        dir.create(paste0(dir_div,"/agg-",area),recursive=TRUE,showWarnings=FALSE)
        currentdiv.c<-aggregate(alphadiv, fact=dim, fun=stat, na.rm=TRUE)

        writeRaster(currentdiv.c,file=paste0(dir_div,"/agg-",area,"/maize_alpha_diversity_",stat,"_",area,"km2_",p,".grd"),overwrite=TRUE)
      } # stat in c("mean","sum")
  } # for dim in  c(5,7,10)
} # for p in projects



#################################################################
# FOR EACH SPECIES GET SUM OF BINARY PREDICTIONS BY AGGREGATING
# TO GET COMMUNITY-LEVEL DATA
#################################################################
### under construction...

           # library(vegan)
           # 
           # 
           # 
           # dir_comm<-paste0(dir_out,"/binaries/comm")
           # dir.create(dir_comm,recursive=TRUE,showWarnings=FALSE)
           # p<-projs[1]
           # 
           # 
           # for (p in projs){
           #   for (sp in sp.n)
           #     sp
           #     p
           #     bin<-raster(paste0(dir_binary,"/binary_tc_70threshold",sp,"_",p,".grd"))
           # 
           #     currentdiv.c<-aggregate(bin, fact=10, fun=sum, na.rm=TRUE#, filename=paste0(dir_range,"/",sp,"bin_sum_25km2",p,".grd")
           #                             )
           #     
           #  
           #     currentdiv.c[is.na(currentdiv.c[])] <- 0 
           #     plot(currentdiv.c)
           #     
           #     fw <- focalWeight(currentdiv.c, 0.2, 'circle')
           #     fw <- ifelse(fw == 0, NA, 1)
           #     
           #         # Neighbourhood richness
           #         richness <- function(x, ...) {
           #           length(unique(na.omit(x)))
           #         }
           #     
           #     
           #     richOut <- focal(currentdiv.c, fw, fun=richness, pad=F)
           #     plot(richOut)
           #         
           #     shannonVegan <- function(x, ...) {
           #           diversity(table(x), index="shannon")
           #     }
           #     
           #         # Neighbourhood Shannon Diversity Index
           #         shannon <- function(x, ...) {
           #           cnts <- table(x)
           #           cnts <- cnts / sum(cnts)
           #           -sum(cnts * log(cnts))
           #         }
           #         
           #     shanVegOut <- focal(currentdiv.c, fw, fun=shannon, pad=F) 
           #     plot()
           #   }
           # }


#################################################################
# GET RANGE CHANGE TO
#################################################################

RangeChangeApply <-function(sp) {
  types=c("binary_750_wmean_avg_","binary_tc_70threshold")
  
  for (type in types){
        #type="binary_tc_70threshold"
    
    
        #################################################################
        # USE BINARY TO CACULATE RANGE SIZE CHANGE BY METRIC
        #################################################################
        preds<-Sys.glob(paste0(dir_binary,"/",type,sp,"_*.grd"))
        currentPred <- raster::raster(preds[1])
        f50Pred  <- raster::raster(preds[2])
        f70Pred <- raster::raster(preds[3])
        
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
        df50chg$species<-sp
        df50chg
        dir.create(paste0(dir_rangedf,"/",type),recursive=TRUE,showWarnings=FALSE)
        dir.create(paste0(dir_rangeR,"/",type),recursive=TRUE,showWarnings=FALSE)
        dir.create(paste0(dir_rangeras,"/",type),recursive=TRUE,showWarnings=FALSE)
        
        write.csv(df50chg,file=paste0(dir_rangedf,"/",type,"/",sp,"_f50_rangesize.csv"))
        writeRaster(rangechange50$Diff.By.Pixel ,file=paste0(dir_rangeras,"/",type,"/",sp,"_2050_diffbypixel.grd"),overwrite=TRUE)
        save(rangechange50,file=paste0(dir_rangeR,"/",type,"/",sp,"_2050_rangechange.Rdata"))
        
        
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
        df70chg$species<-sp
        df70chg
        
        write.csv(df70chg,file=paste0(dir_rangedf,"/",type,"/",sp,"_f70_rangesize.csv"))
        writeRaster(rangechange70$Diff.By.Pixel ,file=paste0(dir_rangeras,"/",type,"/",sp,"_2070_diffbypixel.grd"),overwrite=TRUE)
        save(rangechange70,file=paste0(dir_rangeR,"/",type,"/",sp,"_2070_rangechange.Rdata"))
        
        library(gtable)
        library(grid)
        library(ggplot2)
        
        
        
        ## plot binary differences by pixel for both future periods
        chg<-Sys.glob(paste0(dir_rangeras,"/",sp,"*.grd")) # get ensemble files out
        chgras<-raster(chg[1])
        val <- getValues(chgras)
        xy <- as.data.frame(xyFromCell(chgras,1:ncell(chgras)))
        xy <- cbind(xy,val)
        
        cols<-c("-2"="firebrick", "-1"="royalblue4","0"="cornsilk","1"="skyblue2")
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
        g2 <- ggplotGrob(change)
        g3 <- ggplotGrob(change2)
        g <- rbind(g2, g3, size = "first")
        g$widths <- unit.pmax(g2$widths, g3$widths)
        grid.newpage()
        
        ggsave(paste0(sp,"-",type,"change_by_pixels.png"),width=10, height=10, units="in", plot = grid.draw(g), device = "png",dpi=600, path = paste0(dir_rangeras,"/",type),
               scale = 1)
        
  }
} # for sp in sp.n



library(Rmpi)
library(parallel)
library(snowfall)
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

#ncpus<-ntasks*nnodes
ncpus<-detectCores()
print(paste0(ncpus," CPUs"))

sfInit( parallel=TRUE, cpus=ncpus, type="MPI",
        slaveOutfile = paste0(Sys.getenv("SLURM_JOB_ID"),"_sfInit.log"))

# # Export packages to snowfall
print("Exporting packages to cluster")
sfLibrary('biomod2', character.only = TRUE)
sfLibrary('raster', character.only = TRUE)
sfLibrary('gtable', character.only = TRUE)
sfLibrary('grid', character.only = TRUE)
sfLibrary('ggplot2', character.only = TRUE)
sfLibrary('rasterVis', character.only = TRUE)
sfLibrary('sp', character.only = TRUE)

# # Export environment to snowfall
sfExportAll()
sfExport("RangeChangeApply")

sfLapply(sp.n, RangeChangeApply)

sfStop()



#################################################################
# GET ALPHA DIVERSITY
#################################################################

# # following https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/biomod2/inst/doc/Multi_species_computation.pdf?revision=315&root=biomod
#
# dir_binary<-paste0(dir_out,"/binaries")
# for (p in projects){
#   dir_range<-paste0(dir_out,"/range-size")
#   
#   dir_binary<-paste0(dir_out,"/binaries")
#   
#   bin_grds<-Sys.glob(paste0(dir_binary,"/*",p,"*.grd")) # get ensemble files out
#   bin_stack<-stack(bin_grds)
#   
#   alpha_div<-raster::stackApply(bin_stack,indices=c(rep(1,raster::nlayers(bin_stack))),fun=sum)
#   writeRaster(alpha_div ,file=paste0(dir_range,"/",p,"_alpha_diversity.grd"),overwrite=TRUE)
# }


projs<-c("current","rcp85_50","rcp85_70")


