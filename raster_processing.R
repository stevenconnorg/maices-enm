library(raster)
library(XML)
require(rgdal)
# define directories

root<-"E:/thesis"
setwd(root)

dir_dat<-paste0(root,"/01_data")
dir_R<-paste0(root,"/02_R")
dir_out<-paste0(root,"/03_output")
dir_figs<-paste0(root,"/04_figs")
dir_lit<-paste0(root,"/05_lit")
dir_comp<-paste0(root,"/06_comp")
dir_presentations<-paste0(root,"/07_pres")

dir_maices<-paste0(dir_dat,"/maices")
dir_ind<-paste0(dir_dat,"/ind")
dir_topo<-paste0(dir_dat,"/topo")
# 
dir_clim<-paste0(dir_dat,"/clim")
dir_pres<-paste0(dir_clim,"/present")
dir_fut<-paste0(dir_clim,"/future")

dir_p.mosaics<-paste0(dir_pres,"/2.0/")
dir_f.mosaics<-paste0(dir_fut,"/1.4/")

dir_land<-paste0(dir_dat,"/land-cover")

dir_stacks<-paste0(dir_dat,"/stacks/")
dir_envirem<-paste0(dir_dat,"/envirem/")

folders<-as.list(ls())
i<-folders[[3]]
# tk work on function to create these directories
for (i in 1:length(folders))  { 
  f<-folders[[i]]
  folder<-get(f)
  dir.create(folder,recursive=TRUE) 
} 
dir.create(dir_stacks)

# load(paste0(dir_clim,"/raster_processing.RData"))

ls()
# 


# extent and resolution for final rasters:
# class       : Extent 
# xmin        : -117.625 
# xmax        : -86.20833 
# ymin        : 14.025 
# ymax        : 33.225 
# nrow=2304,ncol=3770

# template raster 
r<-raster(nrow=2304,ncol=3770)
extent(r)<-c(-117.625,-86.20833,14.025,33.225)
r


##############################################
# get soil data from FAO harmonilzed world soil db v1.2
##############################################

urlsoil<-c("http://www.fao.org/soils-portal/soil-survey/soil-maps-and-databases/harmonized-world-soil-database-v12/en/")
docsoil<-htmlParse(urlsoil)
#get <a> nodes.
Anodes<-getNodeSet(docsoil,"//a")
grep("*.asc*",Anodes,value=T)
#make the full url
urls<-grep("asc",sapply(Anodes, function(Anode) xmlGetAttr(Anode,"href")),value=TRUE)
urls<-paste0("http://www.fao.org/",urls)
filename<-basename(urls)
library(tools)
ext<-file_ext(filename)
asc<-paste0(".",ext[1])
file<-gsub(asc,"",filename)
landfiles<-paste0(dir_land,"/",filename)
mapply(function(x,y) download.file(x,y),urls,files)

r<-raster(nrow=2304,ncol=3770)
extent(r)<-c(-117.625,-86.20833,14.025,33.225)


bbox<-bbox(r)
for(i in landfiles){ ##140:173
  name<-gsub(".asc","",basename(i))
  ras<-raster(i)
  ri<-resample(ras,r)
  ri<-round(ri,0)
  croplay<-crop(ri,bbox) ##filtered
  writeRaster(croplay,filename=paste0(dir_land,"/crop/crop_",name,".grd"),overwrite=TRUE)
  do.call(file.remove,list(list.files(pattern="temp*"))) 
}

# rename and reorder stack layers

landfiles<-list.files(paste0(dir_land,"/crop"),full.names = T,pattern = ".grd")
landstack<-stack(landfiles)
names(landstack)
names(landstack)<-c("TotalCult","IrrCult","Rain-fedCult","Forested","Grass/Woodland","Barren","SQ1NutAvail","SQ2NutRetn","SQ3RootingCond","SQ4OxyAvailRoots","SQ5ExcessSalts","SQ6Toxicity","SQ7Workability","Urban","Water")
#plot(landstack)
landcoverstack<-landstack[[c(1:6,14:15)]] #extract land cover layers
writeRaster(landcoverstack,filename=paste0(dir_stacks,"/FAOlandstack.grd"),bylayer=FALSE,format="raster",overwrite=TRUE)

soilqualstack<-landstack[[7:13]] # extract soil quality variables
library(raster)
for (i in 1:nlayers(soilqualstack)){
  soilqualstack[[i]]<-raster::as.factor(soilqualstack[[i]]) # assign soil quality variables as factors
  
}
writeRaster(soilqualstack,filename=paste0(dir_stacks,"/FAOsoilstack.grd"),bylayer=FALSE,format="raster",overwrite=TRUE)

##############################################
# Download WorldClim altitude Rasters by tile
##############################################
tiles<-(c("11","12","13","21","22","23"))
for (i in tiles){
  download.file(paste0("http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/alt_",i,"_tif.zip"),destfile = paste0(dir_topo,"/alt_",i,".zip"),method = "wget")
  zipF<-paste0(dir_topo,"/alt_",i,".zip") # lets you choose a file and save its file path in R (at least for windows)
  outDir<-paste0(dir_topo,"/alt_",i) # Define the folder where the zip file should be unzipped to 
  unzip(zipF,exdir=outDir)  # unzip your file 
}

# altrasters<-list.files(dir_topo,pattern=".bil",full.names = T,recursive=T)

altrasters<-lapply(list.files(dir_topo,pattern=".tif",full.names = T,recursive=T),raster)
altrasters$fun<-mean
alt.mosaic<-do.call(mosaic,altrasters)
crs(alt.mosaic)<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
names(alt.mosaic)<-"altitude"
alt.crop<-crop(alt.mosaic,r)
writeRaster(alt.crop,filename=paste0(dir_topo,"/alt_cropped.grd"),format="raster",overwrite=T)

alt.crop<-raster(paste0(dir_topo,"/alt_cropped.grd"))
names(alt.crop)<- "elev"
terrain<-terrain(alt.crop, opt=c('aspect','slope','Roughness'), unit='degrees', neighbors=8)
names(terrain)
writeRaster(terrain, paste0(dir_topo,"/crop/",c('aspect','slope','Roughness'),".grd"), bylayer=TRUE, format='raster', overwrite=T)
topostack<-stack(terrain,alt.crop)
writeRaster(topostack, paste0(dir_stacks,"topostack.grd"), bylayer=FALSE, format='raster', overwrite=T)


##############################################
# download and unzip present Wordclim 2.0 climate data
##############################################
setwd(dir_p.mosaics)
url<-c("http://worldclim.org/version2")
doc<-htmlParse(url)
#get <a> nodes.
Anodes<-getNodeSet(doc,"//a")

#make the full url
urls<-grep("30s",sapply(Anodes, function(Anode) xmlGetAttr(Anode,"href")),value=TRUE)

#Select the files of interest 
r<-regexpr("wc.*?zip",urls)
files<-regmatches(urls,r)
mapply(function(x,y) download.file(x,y),urls,files)
lapply(grep(".zip",files, value=TRUE),unzip)

# alternatively with a cluster
# get inside present clim data mosaic directory
# setwd(dir_p.mosaics)

# library(parallel)
# cl <- parallel::makeCluster(detectCores())
# clusterMap(cl, download.file, url = urls, destfile = paste0(dir_p.mosaics,"/",files))

# stopCluster(cl)

#####
# Present clim data Raster Stack Cropping
dir.create(paste0(dir_p.mosaics,"/crop"))

# read in present data, stack and write to RData file
pres_rasts<-list.files(paste0(dir_p.mosaics),pattern="\\.tif$", full.names=TRUE)

presfullstack<-stack(pres_rasts)

# write to file
save(presfullstack,file=paste0(dir_stacks,"/present_fullstack.RData"))
writeRaster(presfullstack, paste0(dir_stacks,"/present_fullstack.grd"), bylayer=FALSE, format='GTiff',overwrite=FALSE)

dir.create(paste0(dir_p.mosaics),"/crop")

bbox<-extent(r)
# i<-pres_rasts[1]
pres_rasts<-list.files(paste0(dir_p.mosaics),pattern="srad", full.names=TRUE)

for(i in pres_rasts){ ##140:173
  name<-gsub(".tif","",basename(i))
  ras<-raster(i)
  croplay<-crop(ras,bbox) ##filtered
  writeRaster(croplay,filename=paste0(dir_p.mosaics,"/crop/crop_",name,".tif"),overwrite=TRUE)
  do.call(file.remove,list(list.files(pattern="temp*"))) 
}

save.image(paste0(dir_R,"/raster_processing.RData"))

#  make present crop raster stack
pres_croprasts<-list.files(paste0(dir_p.mosaics,"/crop"),pattern="\\.tif$", full.names=TRUE)

prescropstack<-stack(pres_croprasts)

# write to file
save(prescropstack,file=paste0(dir_stacks,"/present_cropstack.RData"))
writeRaster(prescropstack, paste0(dir_stacks,"/present_cropstack.grd"), bylayer=FALSE, format='raster',overwrite=T)



#### DOWNLOAD FUTURE DATA FOR RCP 85 ####
# from https://github.com/ClimateActionUCI/datasets/blob/master/get_CMIP5.R

library(dplyr)
library(raster)
AR5<-'AR5 temperature increase projections'
ar5.df<-data.frame('Scenario' = c('RCP2.6','RCP4.5','RCP6.0','RCP8.5'),
                   '2046 - 2065' = c('1.0 (0.4 to 1.6)','1.4 (0.9 to 2.0)','1.3 (0.8 to 1.8)','2.0 (1.4 to 2.6)'),
                   '2081 - 2100' = c('1.0 (0.3 to 1.7)','1.8 (1.1 to 2.6)','1.2 (1.4 to 3.1)','3.7 (2.6 to 4.8)'))
mods<-expand.grid(var=c("tn","tx","pr","bi"),   #tn, tx, pr, or bi, no bi?
                  rcp=c(85), ##26, 45, 60, 85   # rcp
                  model=
                    
                    
                    # following Condo et al. 2011 "Regional climate change scenarios for MÃÂ©xico." AtmÃÂ³sfera 24(1), 125-140 (2011)
                    c("CC", # CCSM4 (Community Climate System Model, UCAR)- new version of CCSM-30 
                      "MC", # MIROC5 (Model for Interdisciplinary Research on Climate)- new version of MIROC-HI; "Although its good performance and high resolution, MIROC32-HIRES model has an inconvenience: its sensibility is 5.6 ÃÂºC, way higher than the 3 ÃÂºC marked as Ã¢ÂÂbest estimateÃ¢ÂÂ in IPCCÃÂ´s AR4 (Wigley, 2008). " (Conde et al. 2011)
                      "MP", # MPI-ESM-LR (Max-Plank Institute) - per 5th National Communication of Mexico for the United Nations Framework Convention on Climate Change http://atlasclimatico.unam.mx/atlas/Docs/f_escenarios.html
                      "HE", # HADGEM2-ES (Met Office Hadley)per 5th  removed because already downloaded
                      "GF") # GFDL-CM3 (Geophysical Fluid Dynamics Laboratory )
                  ,
                  
                  year=c(50,70), ##50 or 70     # period 2050 or 2070
                  res="30s")                # resolution
  
  mutate(filename=paste0(tolower(model),rcp,var,year),
         url=paste0("http://biogeo.ucdavis.edu/data/climate/cmip5/",res,"/",filename,".zip"))

# set path of where to download future data zips to
outpath<-dir_f.mosaics

# set download.file timeout 
options(timeout=10000)
# make function to download files
dwnldfxn<-function(aurl,filename){
  try(raster:::.download(aurl,filename))
}

# extract urls from mods df
urls<-mods$url

# get zipfiles
zipfile<-paste0(outpath,substr(urls,nchar(urls)-12+1,nchar(urls)))


# download each zipfile for each url
mapply(dwnldfxn,aurl=urls,filename=zipfile)

###########################################################
## unzip files and crop
###########################################################

# for cropping, read in some data for mexico to get outline
bbox<-extent(r)

# subset
path = paste0(dir_f.mosaics)
zfs<-list.files(path,pattern="zip",full.names=T)
str(zfs)

dir.create(paste0(dir_f.mosaics,"/crop/ensemble/50"),recursive=TRUE)
dir.create(paste0(dir_f.mosaics,"/crop/ensemble/70"),recursive=TRUE)

# unzip zip directories for each climatology into respective directories, crop, and write layers to 'crop' directory
for(i in zfs[]){ ##140:173
  exdir= gsub(".zip","",i)
  unzip(i,exdir=exdir)  # unzip file
  apatt<-substr(i,nchar(i)-12+1,nchar(i)-4)
  gtifs<-list.files(exdir,pattern=".tif",full.names=T)# [c(2, 6:13, 3:5)]##reorder
  t<-gtifs[1]
  for (t in gtifs){
    gtifras<-raster(t)
    cropped<-crop(gtifras,bbox)
    writeRaster(cropped,format="raster",filename=paste0(dir_f.mosaics,"crop/",gsub(".tif","",basename(t)),".grd"),overwrite=T)
    # unlink(t)
    do.call(file.remove,list(list.files(pattern="temp*"))) 
  }
  do.call(file.remove,list(list.files(pattern="temp*"))) 
}




# ensemble GCMs by monthly means
# get list of all tifs in crop directory, mean ensemble by variable, and write layer to stack
cmip5files<-list.files(paste0(dir_f.mosaics,"/crop"),pattern="*.grd",recursive = FALSE)
for (i in cmip5files){
  clims<-gsub(".grd","",basename(i))
  filter<-paste0(substr(clims,3,nchar(clims)),".grd")
  period<-substr(clims,7,8)
  name <-gsub(".grd","",basename(filter))
  varfiles1<-list.files(paste0(dir_f.mosaics,"/crop"),pattern= filter,recursive = F)
  varfiles2<-paste0(dir_f.mosaics,"/crop/",varfiles1)
  stack<-stack(varfiles2)
  stack2<-stackApply(stack,indices=rep(1:1,length(names(stack))),fun="mean")
  names(stack2)<-paste0(name,"_ensemble")
  writeRaster(stack2,filename=paste0(dir_f.mosaics,"crop/ensemble/",period,"/",name,"_ensemble.grd"),overwrite=TRUE)
  # topofiles<-list.files(paste0(dir_dat,"/topo"),pattern="crop",recursive=TRUE)
  # for (t in topofiles){
  #   file.copy(paste0(dir_dat,"/topo/",t),paste0(dir_f.mosaics,"/crop/ensemble/",period,"/",basename(t)))
  # }
}


#  make  crop raster stack
f50_croprasts<-list.files(paste0(dir_f.mosaics,"crop/ensemble/50/"),pattern="\\.grd$", full.names=TRUE)
f70_croprasts<-list.files(paste0(dir_f.mosaics,"/crop/ensemble/70/"),pattern="\\.grd$", full.names=TRUE)

f50cropstack<-stack(f50_croprasts)
f70cropstack<-stack(f70_croprasts)

# write to file

save(f50cropstack,file=paste0(dir_stacks,"/f50cropstack.RData"))
save(f70cropstack,file=paste0(dir_stacks,"/f70cropstack.RData"))
hdr(f50cropstack, format = "ENVI") # to preserve layer names in other programs with .grd/.gri file types -- uncompressed, so they take a while
hdr(f70cropstack, format = "ENVI")
writeRaster(f50cropstack, paste0(dir_stacks,"/f50cropstack.grd"), bylayer=FALSE, format='raster', overwrite=T)
writeRaster(f70cropstack, paste0(dir_stacks,"/f70cropstack.grd"), bylayer=FALSE, format='raster', overwrite=T)
save.image(paste0(dir_clim,"/raster_processing.RData"))



### READ IN ALL RASTER STACKS, OVERWRITE WITH MATCHING LAYER INDICES
library(raster)
library(quickPlot)
library(usdm)
# read in raster stacks
grds<-list.files(path=dir_stacks,pattern=".grd",full.names = T)

pres_cropstack<-stack(grds[7])
f50cropstack<-stack(grds[3])
f70cropstack<-stack(grds[6])

# how many layers in each stack?
nlayers(pres_cropstack)
nlayers(f50cropstack)
nlayers(f70cropstack)

# get layer names
layerNames(pres_cropstack)
layerNames(f50cropstack)
layerNames(f70cropstack)

# rearrange to have matching layer vars
# bioclimatic, prec
# get preslayer crop stack to match future crop stacks

preslayervecs<-c(73:91,1:12,37:48,25:36)
layerNames(pres_cropstack[[preslayervecs]])
nlayers(pres_cropstack[[preslayervecs]])


pres_cropstack<-pres_cropstack[[preslayervecs]]
layerNames(pres_cropstack)

# make sure raster stack names are the same, formatting first
library(quickPlot)
names(pres_cropstack)<-gsub("crop_wc2.0","",layerNames(pres_cropstack)) # remove prefix
names(pres_cropstack)<-gsub("30s_","",layerNames(pres_cropstack))        # remove suffix
names(pres_cropstack)<-gsub("X_","",layerNames(pres_cropstack))        # remove suffix
layerNames(pres_cropstack)

names(f50cropstack)
names(f50cropstack)<-layerNames(pres_cropstack)  # apply presmodstack layer names to future stacks
names(f70cropstack)<-layerNames(pres_cropstack)

# save raster stacks
f50brick<-brick(f50cropstack)
f70brick<-brick(f70cropstack)
presbrick<-brick(pres_cropstack)

hdr(f50brick, format = "ENVI") # to preserve layer names in other programs with .grd/.gri file types -- uncompressed, so they take a while
hdr(f70brick, format = "ENVI")
hdr(presbrick, format = "ENVI")


writeRaster(presbrick, paste0(dir_stacks,"/present_cropstack.grd"), bylayer=FALSE, format='raster', overwrite=T)
writeRaster(f50brick, paste0(dir_stacks,"/f50cropstack.grd"), bylayer=FALSE, format='raster', overwrite=T)
writeRaster(f70brick, paste0(dir_stacks,"/f70cropstack.grd"), bylayer=FALSE, format='raster', overwrite=T)
save.image(paste0(dir_clim,"/raster_processing.RData"))


### MAKE SOME PLOTS
library(raster)
# read in raster stacks
grds<-list.files(path=dir_stacks,pattern=".grd",full.names = T)

pres_cropstack<-stack(grds[7])
f50cropstack<-stack(grds[3])
f70cropstack<-stack(grds[6])

# write rasters to file 
pres_cropstack@title<-"1970-2000"
f50cropstack@title<-"2041-2060"
f70cropstack@title<-"2061-2080"

cropstacks<-list(pres_cropstack,f50cropstack,f70cropstack)
library(rasterVis)
library(quickPlot)

i[1:19]
layerRanges<-list(1:19,20:31,32:43,44:55)
l<-layerRanges[[1]]
i<-cropstacks[[1]]
for (i in cropstacks[i]){
  for (l in layerRanges[[l]]){
    as.integer(l)
    layers<-as.character(l)
    x<-i@title
    label<-strsplit(layerNames(i[[l]]), "[_]")[[1]][1]
    png(filename=paste0(dir_figs,"/",x,"_",label,".png"))
    plot(i[[l]]) # prcp
    dev.off()
  }  
}


for (i in cropstacks){
  x<-i@title
  biolabel<-strsplit(layerNames(i[[1:19]]), "[_]")[[1]][1]
  preclabel<-strsplit(layerNames(i[[20:31]]), "[_]")[[1]][1]
  tminlabel<-strsplit(layerNames(i[[32:43]]), "[_]")[[1]][1]
  tmaxlabel<-strsplit(layerNames(i[[44:55]]), "[_]")[[1]][1]

  png(filename=paste0(dir_figs,"/",x,"_",biolabel,".png"))
  par(mfrow=c(4,5))
  plot(i[[1:19]])
  dev.off()
  
  png(filename=paste0(dir_figs,"/",x,"_",preclabel,".png"))
  par(mfrow=c(3,4))
  plot(i[[20:31]]/10) # prcp
  dev.off()
  
  png(filename=paste0(dir_figs,"/",x,"_",tminlabel,".png"))
  par(mfrow=c(3,4))
  plot(i[[32:43]]/10) # tmin
  dev.off()
  
  png(filename=paste0(dir_figs,"/",x,"_",tmaxlabel,".png"))
  par(mfrow=c(3,4))
  plot(i[[44:55]]/10) # tmax
  dev.off()
}

f50cropstack@title<-"2041-2060"
f70cropstack@title<-"2061-2080"
fcropstacks<-list(f50cropstack,f70cropstack)
#fcropstacks
for (cropstack in fcropstacks){
  preclabel<-strsplit(layerNames(cropstack[[20:31]]), "[_]")[[1]][1]
  tminlabel<-strsplit(layerNames(cropstack[[32:43]]), "[_]")[[1]][1]
  tmaxlabel<-strsplit(layerNames(cropstack[[44:55]]), "[_]")[[1]][1]

  # get max temp diff
  tmaxdiff<-((f70cropstack[[44:55]]/10)-((pres_cropstack[[44:55]])))
  png(filename=paste0(dir_figs,"/",tmaxlabel," (C) 1970 - 2000 avg. to ",cropstack@title," avg. Difference, RCP 8.5.png"))
  rasterVis::levelplot(tmaxdiff,main=paste0(tmaxlabel," (C) 1970 - 2000 avg. to ",cropstack@title," avg. Difference, RCP 8.5"))
  dev.off()
  
  # get max temp diff
  tmindiff<-(cropstack[[32:43]]/10)-(pres_cropstack[[32:43]])
  png(filename=paste0(dir_figs,"/",tminlabel," (C) 1970 - 2000 avg. to ",cropstack@title," avg. Difference, RCP 8.5.png"))
  rasterVis::levelplot(tmindiff,main=paste0(tminlabel," (C) 1970 - 2000 avg. to ",cropstack@title," avg. Difference, RCP 8.5"))
  dev.off()
  
  # get prcp diff
  pcpdiff<-(cropstack[[20:31]])-(pres_cropstack[[20:31]])
  png(filename=paste0(dir_figs,"/",pclabel," (mm) 1970 - 2000 avg. to ",cropstack@title," avg. Difference, RCP 8.5.png"))
  rasterVis::levelplot(pcpdiff,main=paste0(pclabel," 1970 - 2000 avg. to ",cropstack@title," avg. Difference, RCP 8.5"))
  dev.off()
  
}

#####################################
# Generate envirem data
# write rasters to file 
library(raster)
grds<-list.files(path=paste0(dir_pres,"/2.0/crop"),pattern=".tif",full.names = T)

presbio_cropstack<-stack(grep(grds,pattern = "bio",value = T))
presprec_cropstack<-stack(grep(grds,pattern = "prec",value = T))
prestmin_cropstack<-stack(grep(grds,pattern = "tmin",value = T))
prestmax_cropstack<-stack(grep(grds,pattern = "tmax",value = T))

# WorldClim 2.0 has degrees C as temp units
# Worldclim 1.4 has degrees C * 10 as temp units
# make them all have same unit across stacks
<<<<<<< HEAD

prestmin_cropstack<-prestmin_cropstack*10
prestmax_cropstack<-prestmax_cropstack*10

pres_cropstack<-stack(presbio_cropstack,presprec_cropstack,prestmin_cropstack,prestmax_cropstack)
#plot(pres_cropstack[[20:31]])
#plot(f50cropstack[[20:31]])
#plot(f70cropstack[[20:31]])

=======

prestmin_cropstack<-prestmin_cropstack*10
prestmax_cropstack<-prestmax_cropstack*10

pres_cropstack<-stack(presbio_cropstack,presprec_cropstack,prestmin_cropstack,prestmax_cropstack)
#plot(pres_cropstack[[20:31]])
#plot(f50cropstack[[20:31]])
#plot(f70cropstack[[20:31]])

>>>>>>> 7e624be3be034466e77ee95b6955d434d2cde684
plot(pres_cropstack[[20:31]])
plot(f50cropstack[[20:31]])
plot(f70cropstack[[20:31]])

library(quickPlot)
names(pres_cropstack)<-gsub("crop_wc2.0","",layerNames(pres_cropstack)) # remove prefix
names(pres_cropstack)<-gsub("30s_","",layerNames(pres_cropstack))        # remove suffix
names(pres_cropstack)<-gsub("X_","",layerNames(pres_cropstack))        # remove suffix
layerNames(pres_cropstack)

pres_cropstack@title<-"1970-2000"
f50cropstack@title<-"2041-2060"
f70cropstack@title<-"2061-2080"

cropstacks<-list(pres_cropstack,f50cropstack,f70cropstack)

library(envirem)
library(quickPlot)
library(stringr)
for (i in list(pres_cropstack,f50cropstack,f70cropstack)){
  bio<-i[[(grep(layerNames(pres_cropstack),pattern = "bio",value = T))]]
  tmin<-i[[(grep(layerNames(pres_cropstack),pattern = "tmin",value = T))]]
  tmax<-i[[(grep(layerNames(pres_cropstack),pattern = "tmax",value = T))]]
  prec<-i[[(grep(layerNames(pres_cropstack),pattern = "prec",value = T))]]
  #tmin <- i[[32:43]] 
  #tmax <- i[[44:55]] 
  #prec <- i[[20:31]]
  #bio <- i[[1:19]]
  tmean<-mosaic(tmin,tmax,fun=mean)
  names(tmean)<-paste0("tmean_",str_pad(rep(01:12), 2, pad = "0"))
  tmean

  sradFiles<-list.files(paste0(dir_pres,"/2.0/crop"),pattern="srad",full.names = T)
  srad<-stack(sradFiles)
  names(srad)<-paste0("srad_",str_pad(rep(01:12), 2, pad = "0"))
  

  p<-stack(bio,prec,srad,tmax,tmean,tmin)
  p@title<-i@title
  names(p)<-c(paste0("bio_",rep(01:19)),
              paste0("prec_",rep(01:12)),
              paste0("et_solrad_",rep(01:12)),
              paste0("tmax_",rep(01:12)),
              paste0("tmean_",rep(01:12)),
              paste0("tmin_",rep(01:12)))
  
  
  dir<-paste0(dir_stacks,i@title)
  dir.create(dir,showWarnings = FALSE)
  writeRaster(p, paste0(dir,"/",names(p),".tif"), bylayer=TRUE, format='GTiff', overwrite=T)

  generateRasters(var='all', maindir=dir, resName=NULL, timeName=p@title, outputDir=paste0(dir_dat,"/envirem"),
                  rasterExt = ".tif", nTiles = 1, overwriteResults = TRUE,
                  outputFormat = "GTiff", tempDir = "~/temp", gdalinfoPath = NULL,
                  gdal_translatePath = NULL)
}

 env<-stack(list.files(path="E:\\thesis\\01_data\\envirem\\2061-2080\\",pattern=".tif",full.names=T))
plot(env)
#biolabel i[[1:19]]
#prec i[[20:31]]
#tmin i[[32:43]]
#tmax i[[44:55]]

###########################################################
###       MODELLING STACK VARIABLE SELECTION            ###
###########################################################
library(raster)
# read in new cropstacks with all calculated variables

for (cropstack in cropstacks){
  print(cropstack@title)
  biostack<-raster::stack(list.files(path=paste0(dir_stacks,cropstack@title,"/"),pattern="bio",full.names = T))
  enviremstack<-raster::stack(list.files(path=paste0(dir_envirem,cropstack@title,"/"),pattern=".tif",full.names = T))
  
  ### CLEAN UP LAYER NAMES AFTER ENVIREM ###
  # get biostack names back to equal lengths
  names(biostack)[nchar(names(biostack))==5] <- sub("(.{4})(.*)", "\\10\\2", layerNames(biostack)[nchar(layerNames(biostack))==5])
  names(biostack)
  # and reorder
  biostack<-biostack[[c(1,12:19,2:11)]]
  
  names(biostack)
  names(enviremstack)
  
  full_stack<-stack(biostack,enviremstack)
  
  names(full_stack)
  # change envirem time id to match cropstack title with hyphen
  time<-gsub("-",".",cropstack@title)
  # remove id's added by envirem
  names(full_stack)<-gsub(paste0("X",time,"__")
                          ,""
                          ,names(full_stack)
                          )
  
  # set 19 bioclim vars with real name
  names(full_stack)     <-c("BIO1AnnMeanTemp",
                            "BIO2MeanDiurnalRange",
                            "BIO3Isotherm",
                            "BIO4TSeas",
                            "BIO5TWarmestMonth",
                            "BIO6MinTColdestMonth",
                            "BIO7TAnnRange",
                            "BIO8MeanTWettestQ",
                            "BIO9MeanTDriestQ",
                            "BIO10MeanTWarmestQ",
                            "BIO11MeanTColdestQ",
                            "BIO12AnnPrec",
                            "BIO13PrecWettestMonth",
                            "BIO14PrecDriestMonth",
                            "BIO15PrecSeas-COV",
                            "BIO16PrecWettestQ",
                            "BIO17PrecDriestQ",
                            "BIO18PrecWarmestQ",
                            "BIO19PrecColdestQ",
                            "EVMannPET",
                            "EVMthornthwaiteAI",
                            "EVMclimaticMI",
                            "EVMcontinentality",
                            "EVMembergerQ",
                            "EVMgrowingDegDays0",
                            "EVMgrowingDegDays5",
                            "EVMmaxTempColdest",
                            "EVMminTempWarmest",          
                            "EVMmonthCountByT10",
                            "EVMPETColdestQ",
                            "EVMPETDriestQ",
                            "EVMPETseas",
                            "EVMPETWarmestQ",
                            "EVMPETWettestQ",
                            "EVMthermicityIndex"         
  )
  
  writeRaster(full_stack, paste0(dir_stacks,"/",cropstack@title,"_biostack",".grd"), bylayer=FALSE, format='raster', overwrite=T)
  
}



pres_biostack<-stack(paste0(dir_stacks,"/1970-2000_biostack.grd"))
landstack<-stack(paste0(dir_stacks,"/FAOlandstack.grd"))
topostack<-stack(paste0(dir_stacks,"/topostack.grd"))

fullpresstack<-stack(pres_biostack,landstack,topostack)

# create correlation matrix of bioclimatic variables
library(corrplot)
fullpresstack[is.na(fullpresstack)] <- 0

mat<-as.matrix(fullpresstack)
colnames(mat)<-names(fullpresstack)

#View(mat)
cormat<-cor(mat)
# visualize correlation matrix

png(filename=paste0(dir_figs,"/biovars_cormat.png"),
    width = 1000, height = 1000)
par(mfrow=c(1,1))
corrplot(cormat)
dev.off()



write.csv(mat,file=paste0(dir_out,"/biovars_matrix.csv"))
#mat<-read.csv(file=paste0(dir_out,"/biovars_matrix.csv"))
write.csv(cormat,file=paste0(dir_out,"/biovars_corr_matrix.csv"))


library(usdm)
# get variance inflation factors of pres_cropstack vars
presvif<-usdm::vif(fullpresstack)

layers<-c() # remove layers with high vif, selecting variables appropriate to species (e.g.: maize)
presviflay<-mat[1:50]


presvifstep<-vifstep(fullpresstack,th=10)

presvifstep<-vifstep(pres_biostack,th=10)

presvifcor<-vifcor(fullpresstack,th=0.75)
save.image(file=paste0(dir_out,"/vif.RData"))



# remove multicollinearity of full stack3
#library(SpaDES)
#library(virtualspecies)
#coll_vars_pres<-virtualspecies::removeCollinearity(pres_cropstack,multicollinearity.cutoff = 0.45, select.variables = FALSE, sample.points = F,plot = TRUE)
#save.image(paste0(dir_clim,"/raster_processing.RData"))
#print(coll_vars_pres)
#coll_vars_pres


# manually select raster layers to retain from climate pca, including other variables of interest
presmodstack<-raster::stack(paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_01.tif"),
                            paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_02.tif"),
                            paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_03.tif"),
                            paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_05.tif"),
                            paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_06.tif"),
                            paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_15.tif"),
                            paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_18.tif"),
                            paste0(dir_topo,"/alt_cropped.grd"),
                            paste0(dir_ind,"/pob-ind.grd"),
                            paste0(dir_land,"/crop-land-cover.tif")
                            
)

f50modstack<-raster::stack(paste0(dir_f.mosaics,"crop/ensemble/50/85bi501_ensemble.grd"),
                           paste0(dir_f.mosaics,"crop/ensemble/50/85bi502_ensemble.grd"),
                           paste0(dir_f.mosaics,"crop/ensemble/50/85bi503_ensemble.grd"),
                           paste0(dir_f.mosaics,"crop/ensemble/50/85bi505_ensemble.grd"),
                           paste0(dir_f.mosaics,"crop/ensemble/50/85bi506_ensemble.grd"),
                           paste0(dir_f.mosaics,"crop/ensemble/50/85bi5015_ensemble.grd"),
                           paste0(dir_f.mosaics,"crop/ensemble/50/85bi5018_ensemble.grd"),
                           paste0(dir_topo,"/alt_cropped.grd"),
                           paste0(dir_ind,"/pob-ind.grd")
                           
                           
)

f70modstack<-raster::stack(paste0(dir_f.mosaics,"crop/ensemble/70/85bi701_ensemble.grd"),
                           paste0(dir_f.mosaics,"crop/ensemble/70/85bi702_ensemble.grd"),
                           paste0(dir_f.mosaics,"crop/ensemble/70/85bi703_ensemble.grd"),
                           paste0(dir_f.mosaics,"crop/ensemble/70/85bi705_ensemble.grd"),
                           paste0(dir_f.mosaics,"crop/ensemble/70/85bi706_ensemble.grd"),
                           paste0(dir_f.mosaics,"crop/ensemble/70/85bi7015_ensemble.grd"),
                           paste0(dir_f.mosaics,"crop/ensemble/70/85bi7018_ensemble.grd"),
                           paste0(dir_topo,"/alt_cropped.grd"),
                           paste0(dir_ind,"/pob-ind.grd")       
)

       
       


library(raster)
save(presmodstack,file=paste0(dir_stacks,"/present_modstack.RData"))
writeRaster(presmodstack, paste0(dir_stacks,"/present_modstack.grd"), bylayer=FALSE, format='raster',overwrite=TRUE)

save(f50modstack,file=paste0(dir_stacks,"/f50_modstack.RData"))
writeRaster(f50modstack, paste0(dir_stacks,"/f50_modstack.grd"), bylayer=FALSE, format='raster',overwrite=TRUE)

save(f70modstack,file=paste0(dir_stacks,"/f70_modstack.RData"))
writeRaster(f70modstack, paste0(dir_stacks,"/f70_modstack.grd"), bylayer=FALSE, format='raster',overwrite=TRUE)


save.image(file=paste0(dir_clim,"/raster_processing.RData"))
       

# make sure raster stack names are the same, formatting first
library(quickPlot)
names(presmodstack)<-gsub("crop_wc2.0_","",layerNames(presmodstack)) # remove prefix
names(presmodstack)<-gsub("_30s","",layerNames(presmodstack))        # remove suffix
names(f50modstack)<-layerNames(presmodstack)  # apply presmodstack layer names to future stacks
names(f70modstack)<-layerNames(presmodstack)

title(main = "Average (1970 - 2000) Conditions")
raster::spplot(presmodstack)

modstacks<-c(presmodstack,f50modstack,f70modstack)
for (stack in modstacks){
  for (lay in stack) {
    png(file=paste0(dir_figs,"/",names(lay),"_",paste0(stack)))
    level.plot(layer,main=names(lay))
  }
}
plot(f50modstack)
plot(f70modstack)

