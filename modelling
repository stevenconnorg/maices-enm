getwd()
setwd("E:\\thesis\\02_R")
source(paste0(getwd(),"\\02_R\\functions.R"))

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

# establish directories
root<-"E:\\thesis\\02_R"
setwd(root)
dir_dat<-paste0(root,"/01_data")
dir_R<-paste0(root,"/02_R")
dir_out<-paste0(root,"/03_output")
dir_figs<-paste0(root,"/04_figs")
dir_bm<-paste0(root,"/06_biomod")
dir_bmz<-paste0(root,"/07_biomodez")


# climate data root
dir_clim<-paste0(dir_dat,"/clim")

# present climate data directory
dir_pres<-paste0(dir_clim,"/present")

# future climate data directory
dir_fut<-paste0(dir_clim,"/future")

# present
dir_p.mosaics<-paste0(dir_pres,"/mosaics")
dir_p.raw<-paste0(dir_pres,"/raw")
dir_p.stack<-paste0(dir_pres,"/stack")

# future
dir_f.mosaics<-paste0(dir_fut,"/mosaics")
dir_f.raw<-paste0(dir_fut,"/raw")
dir_f.stack<-paste0(dir_fut,"/stack")


#####################################
#       BIOMOD formatting Data      #
#####################################

# set wd inside biomod subdirectory
setwd(dir_bm)
load(paste0(dir_out,"/stack.R"))
load(paste0(dir_out,"/cropstack.R"))
load(paste0(dir_out,"/cropstack-mcl-rm.R"))
# load(paste0(dir_bm,"/.RData"))

pam<-read.csv(paste0(dir_out,"/pam.csv"))
pa<-data.frame(pam)
pa
maxentjar<-paste0(dir_bm,"/maxent.jar")

sp.n= c("Arrocillo.Amarillo")     #vector of species name(s)

# BioModApply <-function(sp.n) {
  
  myRespName = sp.n
  myResp <- as.numeric(pa[,myRespName])
  myExpl<-presmodstack
  myExplFuture50<-fut50modstack
  myExplFuture70<-fut70modstack
  myRespXY = pa[,c('Longitude.x.','Latitude.y.')]
  
  
  # format input data for biomod
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName
                                       )
  
  # define model options
  
  
  # print default biomod options 
  # default_ModelOptions <-BIOMOD_ModelingOptions()
  # print(default_ModelOptions)
  
  # edit default options accordingly
  BIOMOD_ModelOptions <- BIOMOD_ModelingOptions(GLM = list( type = 'quadratic',       # tk polynomial ?
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
                                                  
                                                 # GAM = list( algo = 'GAM_mgcv',
                                                 #              type = 's_smoother',
                                                #             k = -1,
                                                #             interaction.level = 0,
                                                #              myFormula = NULL,
                                                #             family = binomial(link = 'logit'),
                                                #              method = 'GCV.Cp',
                                                #              optimizer = c('outer','newton'),
                                                #             select = FALSE,
                                                #             knots = NULL,
                                                #             paraPen = NULL,
                                                #             control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07, maxit = 200, trace = FALSE
                                                #                            , mgcv.tol = 1e-07, mgcv.half = 15, rank.tol = 1.49011611938477e-08
                                                #                            , nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
                                                #                            , optim = list(factr=1e+07), newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
                                                #                           , outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE, scale.est = fletcher
                                                #                            , edge.correct = FALSE) ),
                                                  
                                                  
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
  models = c("GLM","GBM","ANN","SRE","CTA","RF","MARS","FDA","MAXENT.Phillips",'MAXENT.Tsuruoka'), 
  models.options = BIOMOD_ModelOptions, 
  NbRunEval=5,
  DataSplit=65,
  Prevalence=NULL,
  VarImport=5,
  models.eval.meth = c( 'KAPPA', 'TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = TRUE,
  modeling.id = paste(myRespName,"FirstModeling",sep=""))
  
maxentBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData, 
  models = c("MAXENT.Phillips"), 
  models.options = BIOMOD_ModelOptions, 
  NbRunEval=5,
  DataSplit=65,
  Prevalence=NULL,
  VarImport=5,
  models.eval.meth = c( 'KAPPA', 'TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = TRUE,
  modeling.id = paste(myRespName,"MaxEntModeling",sep=""))

  myBiomodModelEval <- get_evaluations(myBiomodModelOut)          #GETTING MODEL EVALUATIONS FOR PERFORMANCE TESTING
  print("Done Running Models")

  myBiomodModelEval <- get_evaluations(myBiomodModelOut)          #GETTING MODEL EVALUATIONS FOR PERFORMANCE TESTING
  print("Done Running Models")
  
  plot(get_predictions(myBiomodModelOut))
  
  
  # let's print the ROC scores of all selected models
  # myBiomodModelEval["ROC","Testing.data",,,]
  
  # print variable importances
  # get_variables_importance(myBiomodModelOut)
  
  
  setwd(dir_bm)
  GBMs<- BIOMOD_LoadModels(myBiomodModelOut, models='GBM')
  # get_formal_model(Arrocillo.Amarillo_AllData_Full_GBM) 
  # summary(get_formal_model(Arrocillo.Amarillo_AllData_Full_GBM) )
  
  myRespPlot2D <- response.plot2(models  = GBMs,
                                 Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                                 show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                                 do.bivariate = FALSE,
                                 fixed.var.metric = 'median',
                                 col = c("blue", "red"),
                                 legend = TRUE,
                                 data_species = get_formal_data(myBiomodModelOut,'resp.var'))
  
  
  # BIOMOD_LoadModels(myBiomodModelOut, models='RF')
  # get_formal_model(Arrocillo.Amarillo_AllData_Full_RF) 
  # summary(get_formal_model(Arrocillo.Amarillo_AllData_Full_RF) )
  
  # model projections
  myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExpl,
    proj.name = 'current',
    selected.models = 'all',
    binary.meth = 'none',
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
    em.by='all',
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
