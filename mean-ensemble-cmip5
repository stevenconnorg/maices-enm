cmip5files<-list.files(paste0(dir_f.mosaics,"/crop"),pattern="*.tif",recursive = FALSE)

for (i in cmip5files){
  clims<-gsub(".tif","",basename(i))
  filter<-paste0(substr(clims,3,nchar(clims)),".tif")
  period<-substr(clims,7,8)
  name <-gsub(".tif","",basename(filter))
  varfiles1<-list.files(paste0(dir_f.mosaics,"/crop"),pattern= filter,recursive = TRUE)
  varfiles2<-paste0(dir_f.mosaics,"/crop/",varfiles1)
  stack<-stack(varfiles2)
  stack2<-stackApply(stack,indices=nlayers(stack),fun="mean")
  writeRaster(stack2,filename=paste0(dir_f.mosaics,"/crop/ensemble/",period,"/",name,"_ensemble.tif"),overwrite=TRUE)
  topofiles<-list.files(paste0(dir_dat,"/topo"),pattern="crop",recursive=TRUE)
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  for (t in topofiles){
    file.copy(paste0(dir_dat,"/topo/",t),paste0(dir_f.mosaics,"/crop/ensemble/",period,"/",basename(t)))
  }
}
