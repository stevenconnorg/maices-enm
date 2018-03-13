# Initialize variables

root<-"/gpfs/home/scg67/thesis"
#root<-"E:/thesis"


# directories
dir_dat<-paste0(root,"/01_data")
dir_R<-paste0(root,"/02_R")

dir_out<-paste0(root,"/03_output")
dir_outf<-paste0(root,"/03_output-full")

dir_figs<-paste0(root,"/04_figs")
dir_lit<-paste0(root,"/05_lit")
dir_comp<-paste0(root,"/06_comp")
dir_presentations<-paste0(root,"/07_pres")
dir_maices<-paste0(dir_dat,"/maices")
dir_ind<-paste0(dir_dat,"/ind")

dir_bm<-paste0(dir_R,"/00_biomod")
dir_bmf<-paste0(dir_R,"/00_biomod-full")

dir_topo<-paste0(dir_dat,"/topo")
dir_clim<-paste0(dir_dat,"/clim")
dir_pres<-paste0(dir_clim,"/present")
dir_fut<-paste0(dir_clim,"/future")
dir_p.mosaics<-paste0(dir_pres,"/2.0/")
dir_f.mosaics<-paste0(dir_fut,"/1.4/")
dir_stacks<-paste0(dir_dat,"/stacks/")

folders<-as.list(ls(),pattern="dir")

# tk work on function to create these directories
for (i in folders)  { 
  folder<-get(i)
  dir.create(folder,recursive=TRUE,showWarnings=FALSE) 
  rm(folder)

} 
rm(i)

metrics = c( 'KAPPA', 'TSS', 'ROC')
projects<-c("proj_current","proj_rcp85_50","proj_rcp85_70")
projs<-c("current" , "rcp85_50" ,"rcp85_70")

models = c("MAXENT.Phillips")
