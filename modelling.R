
# establish directories
root<-"E:\\thesis"
setwd(root)

# new directories for biomod
dir_bm<-paste0(dir_R,"/01_biomod")
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
                 "tidyr",
                 "raster",
                 "snowfall"))

Install_And_Load(requiredPackages)

#####################################
#       BIOMOD formatting Data      #
#####################################

# set wd inside biomod subdirectory
setwd(dir_bm)
load(paste0(dir_stacks,"/present_modstack.RData"))
load(paste0(dir_stacks,"/f50_modstack.RData"))
load(paste0(dir_stacks,"/f70_modstack.RData"))
# load(paste0(dir_bm,"/.RData"))

### formating presence absense data, uncomment when on cluster ###
# maices<-readOGR(dsn=paste0(dir_out,"/maices.shp"),layer="maices")
# xy<-cbind(maices@data$Longitud,maices@data$Latitud)
# PA<-letsR::lets.presab.points(xy,maices@data$Raza_prima, resol = 0.01)
# PA<-letsR::lets.presab.points(xy,maices@data$Raza_prima, resol = 0.008333333)
# pam<-PA$Presence_and_Absence_Matrix
# View(pam)
# [pam == 0 ] <- NA
# pa<-data.frame(pam)

# comment out when running cluster
pam<-read.csv(file=paste0(dir_out,"/pam.csv"))


pa<-data.frame(pam)
write.csv(pa,paste0(dir_out,"/pa_dataframe.csv"))

pa
maxentjar<-paste0(dir_bm,"/maxent.jar")

sp.n= c("Tuxpeño"
        #,"Arrocillo.Amarillo"
        )     #vector of species name(s)

library(biomod2)
library(mgcv)
options(max.print=1000000)
library(gbm)
library(dismo)

# BioModApply <-function(sp.n) {
  myRespName = sp.n
  myResp <- as.numeric(pa[,myRespName])
  myExpl<-presmodstack
  myExplFuture50<-f50modstack
  myExplFuture70<-f70modstack
  myRespXY = pa[,c('Longitude.x.','Latitude.y.')]

  
  # format input data for biomod
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       resp.xy = myRespXY,
                                       expl.var = myExpl,
                                       resp.name = myRespName,
                                       PA.nb.rep = 0,
                                       PA.nb.absences = round((sum(!is.na(myResp))/5),0),
                                       PA.strategy = "sre",
                                       PA.sre.quant = 0.95,
                                       na.rm=TRUE
                                       )
  
  # define model options
  
  
  # print default biomod options 
  default_ModelOptions <-BIOMOD_ModelingOptions()
  print(default_ModelOptions)
  
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
                                                              control = list(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001, maxdepth = 25) ),
                                                  
                                                  
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

# modeling
myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData, 
  models = c("GLM","GAM","GBM","ANN","CTA","RF","MARS","FDA","MAXENT.Phillips",'MAXENT.Tsuruoka'), 
  models.options = BIOMOD_ModelOptions, 
  NbRunEval=5,
  DataSplit=60,
  VarImport=5,
  models.eval.meth = c( 'KAPPA', 'TSS', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = TRUE,
  modeling.id = paste(myRespName,"_current"))

print(paste0("Done Running Models for ",sp.n))

          # write data used for modelling
          capture.output(get_formal_data(myBiomodModelOut),
                         file=paste0(dir_out,"/model-data/",myRespName,"_formal_data.txt"))
          ### eval current model
        
          print(paste0("Capturing Model Evaluations for ",sp.n))
          evalmods<-get_evaluations(myBiomodModelOut,as.data.frame=TRUE)
          write.csv(evalmods,file=paste0(dir_out,"/eval/",myRespName,"_formal_models_evaluation.csv"))
          
          ### get variable importance
          modevalimport<-get_variables_importance(myBiomodModelOut,as.data.frame=TRUE)
          write.csv(modevalimport,file=paste0(dir_out,"/var-imp/",myRespName,"_var_imp.csv"))
         
          ### get model summaries
          capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_PA1_Full_ANN",sep="")))))
                         ,file=paste0(dir_out,"/eval/",myRespName,"_ANN_summary.txt"))
          
          capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_PA1_Full_CTA",sep="")))))
                          ,file=paste0(dir_out,"/eval/",myRespName,"_CTA_summary.txt"))
          
          fda1<-get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_PA1_Full_FDA",sep=""))))
          capture.output(fda1$confusion,file==paste0(dir_out,"/eval/",myRespName,"_FDA-conf_summary.txt"))
        
          capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_PA1_Full_GAM",sep="")))))
                         ,file=paste0(dir_out,"/eval/",myRespName,"_GAM_summary.txt"))
          capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_PA1_Full_GBM",sep="")))))
                        ,file=paste0(dir_out,"/eval/",myRespName,"_GBM_summary.txt"))
          
          capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_PA1_Full_GLM",sep="")))))
                         ,file=paste0(dir_out,"/eval/",myRespName,"_GLM_summary.txt"))
        
          capture.output(summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_PA1_Full_MARS",sep="")))))
                        ,file=paste0(dir_out,"/eval/",myRespName,"_MARS_summary.txt"))
          
          # summary(get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_PA1_Full_MAXENT.Phillips",sep="")))))
          maxent_t1<-get_formal_model(get(load(paste(myRespName,"/models/",myRespName,"/",myRespName,"_PA1_Full_MAXENT.Tsuruoka",sep=""))))
        
          myRespPlot2D <- response.plot2(models  = "MAXENT",
                                         Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                                         show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                                         do.bivariate = FALSE,
                                         fixed.var.metric = 'median',
                                         col = c("blue", "red"),
                                         legend = TRUE,
                                         data_species = get_formal_data(myBiomodModelOut,'resp.var'))
          
  # model projections
  myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExpl,
    proj.name = 'current',
    selected.models = 'all',
    binary.meth = c( 'KAPPA', 'TSS','ROC'),
    compress = TRUE,
    clamping.mask = TRUE,
    output.format = '.grd')
  
  
  coutputFolderName <- "proj_current"
  coutputFolder <- paste(myRespName,coutputFolderName,sep="/")      #FINDING DIRECTORY FOR ENSEMBLE MODEL OUTPUT FOR REPORT GENERATION
  print("Done Projecting Models")

  
  mod_proj <- get_predictions(myBiomodProj) 
  plot(mod_proj)
  plot(subset(mod_proj,2), main = "RF projections") # select layer in stack to plot
  
  # ensemble modeling
  myBiomodEM <- BIOMOD_EnsembleModeling(
    modeling.output = myBiomodModelOut,
    chosen.models = 'all',
    em.by='algo',
    eval.metric = c('KAPPA', 'TSS','ROC'),
    eval.metric.quality.threshold = c(0.5,0.5,0.5),
    prob.mean = T,
    prob.cv = T,
    prob.ci = T,
    prob.ci.alpha = 0.05,
    prob.median = T,
    committee.averaging = T,
    prob.mean.weight = T,
    prob.mean.weight.decay = 'proportional' )
  
  get_evaluations(myBiomodEM)
  
  # current ensemble projection
  myBiomodEF <- BIOMOD_EnsembleForecasting(
    EM.output = myBiomodEM,
    projection.output = myBiomodProj)

  myBiomodEF
  EF_stack<- raster::stack(paste0(dir_bm,"/",sp.n,"/proj_current/proj_current_",sp.n,"_ensemble.grd"))
  
  # plot a layer in the ensemble stack
  plot(EF_stack[[6]])
  
  # future projections for rcp 85
  # myBiomodProjFuture70 <- BIOMOD_Projection(
  #   modeling.output = myBiomodModelOut,
  #   new.env = myExplFuture70,
  #   proj.name = 'rcp85_70',
  #   selected.models = 'all',
  #   binary.meth = 'TSS',
  #   compress = 'xz',
  #   clamping.mask = T,
  #   output.format = '.grd')
  
  # future projections for rcp 85
  # myBiomodProjFuture50 <- BIOMOD_Projection(
  #   modeling.output = myBiomodModelOut,
  #   new.env = myExplFuture50,
  #   proj.name = 'rcp85_50',
  #   selected.models = 'all',
  #   binary.meth = 'TSS',
  #   compress = 'xz',
  #   clamping.mask = T,
  #   output.format = '.grd')
  
  
  # f70BiomodEF <- BIOMOD_EnsembleForecasting(
  #   EM.output = myBiomodEM,
  #   projection.output = myBiomodProjFuture70)
  # cat("\n\nExporting Ensemble as ASCII ...\n\n")
  
  # f50BiomodEF <- BIOMOD_EnsembleForecasting(
  #   EM.output = myBiomodEM,
  #   projection.output = myBiomodProjFuture50)
  # cat("\n\nExporting Ensemble as ASCII ...\n\n")
  
  #EXPORTING ENSEMBLE MODEL PROJECTION AS ASCII FOR USE IN OUTSIDE MAPPING SOFTWARE
  gridName = paste(coutputFolderName,myRespName,"ensemble.grd",sep="_")  
  gridDir = paste(myRespName,coutputFolderName,gridName,sep="/")
  MyRaster = raster(paste(myRespName,coutputFolderName,gridName,sep="/"))
  writeRaster(MyRaster,file=paste0(sp.n,"_EnsembleRaster.asc"), format = 'ascii', overwrite = TRUE)
  print("Done Exporting Ensemble as ASCII")

  cat("\n\nExporting Model Plots ...\n\n")
  #EXPORTING PLOTS FOR EACH MODEL & ENSEMBLE TO 'PLOTS' FOLDER IN WD

      model_names <- c("GLM","GBM","GAM","ANN","SRE","CTA","RF","MARS","FDA","MAXENT.Phillips",'MAXENT.Tsuruoka')
  
  # FUNCTION FOR PLOTTING PROJECTIONS FOR EACH MODEL AND PLACEHOLDERS FOR FAILED MODELS
      plot_raster <- function(x,model){
          tryCatch({
              file_name <- paste("plots/", model,".png",sep="")
             png(file_name)
              raster::plot(x, str.grep = model)
        dev.off()
            }, error = function(e){
              file.copy("ParametersAndSettings/no_plot.png",file_name)
              cat("\n\nPlotting for ", model, " failed! Blank file created.\n\n",sep="")
            }
         )
        }
  
      for (i in 1:length(model_names)){
          plot_raster(myBiomodProj,model_names[i])
        }
  
      png(paste0("plots/ensemble_",sp.n,".png"))
      raster::plot(myBiomodEF)
      dev.off()
        print("Done Exporting Model Plots")
  
# }




#####################################################
###############VARIABLE DECISION TREES###############
#####################################################

cat("\n\nBuilding Decision Tree ...\n\n")
p.table <- read.csv(prstbl)                             #TABLE OF PRESENCE/ABSENCE DATA
xcol <- match(xname, names(p.table))                    #FINDING COLUMN NUMBER IN p.table FOR LOGITUDE
ycol <- match(yname, names(p.table))                    #FINDING COLUMN NUMBER IN p.table FOR LATITUDE
pacol <- match(myRespName, names(p.table))              #FINDING COLUMN NUMBER IN p.table FOR SPECIES
xy <- cbind(p.table[xcol], p.table[ycol])               #MAKING NEW OBJECT FOR LONG/LAT ONLY
sp <- SpatialPoints(xy)                                 #CONVERTING TO SPATIALPOINTS OBJECT
xy.extract <- extract(myExpl, sp)                       #EXTRACTING RASTER DATA TO POINTS IN SPATIALPOINTS OBJECT
xy.extract <- cbind(xy.extract, p.table[pacol])         #ADDING PRESENCE/ABSENCE COLUMN TO EXTRACTED RASTER VALUE OBJECT
write.csv(xy.extract, file = "xy_value_extract.csv")    #SAVING EXTRACTED VALUES TO CSV IN WD
xy.extract <- read.csv("xy_value_extract.csv")[,-1]     #IMPORTING CSV BACK IN, FIXING TYPE/CLASS ISSUES
idx_sp <- ncol(xy.extract)                              #FINDING LAST COLUMN WHERE SPECIES PRESENCE/ABSENCE DATA IS LOCATED
dat <- data.frame(xy.extract[,-idx_sp],
                  Species = as.factor(ifelse(xy.extract[,idx_sp] == 0, "absent","present")))        #MAKING NEW SPECIES COLUMN WITH PRESENCE/ABSENCE DATA
xy.tree <- rpart(Species ~ .,dat)                                                                   #BUILDING DECISION TREE FOR PRESENCE/ABSENCE AGAINST ALL EXTRACTED VARIABLE VALUES
dir.create("plots")
png("plots/decision_tree.png")
rpart.plot(xy.tree)            #SAVING PLOT OF DECISION TREE FOR REPORT GENERATION
dev.off()
print("Done Building Decision Tree")

#####################################################
##################GENERATING REPORT##################
#####################################################

cat("\n\nExporting Model Plots ...\n\n")
#EXPORTING PLOTS FOR EACH MODEL & ENSEMBLE TO 'PLOTS' FOLDER IN WD
dir.create("plots")
model_names <- c("GLM","GBM","GAM","ANN","SRE","CTA","RF","MARS","FDA","MAXENT.Phillips",'MAXENT.Tsuruoka')

#FUNCTION FOR PLOTTING PROJECTIONS FOR EACH MODEL AND PLACEHOLDERS FOR FAILED MODELS
plot_raster <- function(x,model){
  tryCatch({
    file_name <- paste("plots/", model,".png",sep="")
    png(file_name)
    raster::plot(x, str.grep = model)
    dev.off()
  }, error = function(e){
    file.copy("ParametersAndSettings/no_plot.png",file_name)
    cat("\n\nPlotting for ", model, " failed! Blank file created.\n\n",sep="")
  }
  )
}

for (i in 1:length(model_names)){
  plot_raster(myBiomodProj,model_names[i])
}

png(paste0("plots/",sp.n,"_ensemble.png"))
raster::plot(myBiomodEF)
dev.off()
print("Done Exporting Model Plots")

#SAVING S4 STRUCTURES FOR REPORT GENERATION
saveRDS(myBiomodModelOut,"myBiomodModelOut.rds")
saveRDS(myBiomodModelEval,"myBiomodModelEval.rds")
variableimportance <- get_variables_importance(myBiomodModelOut)
saveRDS(variableimportance,"variableimportance.rds")
ensembleevaluation <- get_evaluations(myBiomodEM)
saveRDS(ensembleevaluation,"ensembleevaluation.rds")

cat("\n\nLGenerating Report ...\n\n")

#KNITTING RMARKDOWN PDF FROM 'Biomod2_Report.rmd' IN WD
library(rmarkdown)
rmarkdown::render('Biomod2_Report.rmd',output_format=html_document())                           #KNIT HTML REPORT
#Sys.setenv(PATH = paste(Sys.getenv("PATH"),
#    "C:/Program Files/MiKTeX 2.9/miktex/bin/x64", 
#    sep=.Platform$path.sep))                                #ADDING MikTeX to PATH, WILL BE DIFFERENT FOR OTHER LaTeX INSTALLATIONS
#rmarkdown::render('Biomod2_Report.rmd',
#     output_format=pdf_document(latex_engine='xelatex'))    #KNIT PDF REPORT, WILL REQUIRE LaTeX INSTALLATION AND SETUP.
print("Done Generating Report")

################################
# Using Snowfall with lapply over species names
################################

# snowfall initialization
# sfInit(parallel=TRUE, cpus=4)
## Export packages to snowfall
# sfLibrary('biomod2', character.only=TRUE)
## Export variables
# sfExport('myRespXY')
# sfExport('myExpl')
# sfExport('sp.n')
# sfExport('pa')
# sfExport('stack3')
# sfExport('root')
sfExport('BioModApply')
sfExportAll()

# you may also use sfExportAll() to export all your workspace variables
## Do the run
mySFModelsOut <- sfLapply( sp.n, BioModApply)

## stop snowfall
sfStop( nostop=FALSE )



# Sys.glob(file.path(root, sp.n, "current", "R", "*.rdx"))
