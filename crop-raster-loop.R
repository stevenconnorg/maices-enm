# need bbox in environment
# i in list of raster files (pres rasts)

for(i in pres_rasts){ ##140:173
  name<-gsub(".tif","",basename(i))
  ras<-raster(i)
  croplay<-crop(ras,bbox) ##filtered
  writeRaster(croplay,filename=paste0(dir_p.mosaics,"/crop/crop_",name,".tif"),overwrite=TRUE)
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  }
