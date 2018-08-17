# zfs, list of zipfiles

unzipAndcropRasters <- function(zipsfolderPath,removeZip = TRUE,removeUncropped=TRUE)

path = paste0(zipsfolderPath)
zfs<-list.files(path,pattern="zip",full.names=T)
str(zfs)

for(i in zfs){ 
  exdir= gsub(".zip","",i)
  unzip(i,exdir=exdir)  # unzip file
  if (removeZip == TRUE){
	unlink(i)  # remove zip file
	}
  patt<-substr(i,nchar(i)-12+1,nchar(i)-4)
  gtifs<-list.files(exdir,pattern=patt,full.names=T)# [c(2, 6:13, 3:5)]##reorder
  tempstack<-stack(gtifs) ##rasterbrick
  ctempstack<-crop(tempstack,bbox) ##filtered
  writeRaster(ctempstack,bylayer=TRUE,filename=paste0(zipsfolderPath,"/crop/",patt,".tif"))
  if (removeUncropped == TRUE){
	unlink(gtifs)  # remove zip file
	}
  print(paste0("Finished with file ",patt," (",which(zfs==i)," out of ",length(zfs),")"))
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  
  }
