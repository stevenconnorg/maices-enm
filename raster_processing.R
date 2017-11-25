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

dir_stacks<-paste0(dir_clim,"/stacks/")

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

###################################################
# Maize Observation Data Downloading and Cleaning # 
###################################################

download.file("http://www.conabio.gob.mx/informacion/gis/maps/geo/maizngw.zip",destfile = paste0(dir_maize,"/maizngw.zip"),method = "wget")

zipF<-paste0(dir_maices,"maizngw.zip") # lets you choose a file and save its file path in R (at least for windows)
outDir<-paste0(dir_maices) # Define the folder where the zip file should be unzipped to 
unzip(zipF,exdir=outDir)  # unzip your file 


file <- list.files(outDir,pattern = "maizngw",full.names = F) 

#  rename original .shp file name
sapply(files,FUN=function(eachPath){
  file.rename(from=eachPath,to=sub(pattern= paste("maizngw.*"),replacement= "todos-maices.*",eachPath))
})

todos<-readOGR(paste0(dir_dat,"/maices/todos.shp"),layer="todos",use_iconv=TRUE) 

# set encolding to latin 1 for some columns to read special characters
todos@data$Raza_prima<-iconv(todos@data$Raza_prima, from="UTF-8", to="LATIN1")
todos@data$Nom_ent<-iconv(todos@data$Nom_ent, from="UTF-8", to="LATIN1")
todos@data$Nom_mun<-iconv(todos@data$Nom_mun, from="UTF-8", to="LATIN1")

# inspect data
colnames(todos@data) 
str(todos@data)
head(todos@data$Nom_ent)
head(todos@data$Nom_mun)
head(todos@data$NomComun)
length(todos@data[is.na(todos@data$NomComun),])
length(todos@data[is.na(todos@data$Raza_prima),])
unique(todos@data$Raza_prima)


# remove data without prime race information
todos.1<-todos[!is.na(todos@data$Raza_prima),]
# remove data without prime race information
todos.2<-todos.1[todos.1@data$Raza_prima != "ND",]

# remove data where longitude is equal to zero
todos.3<-todos.2[todos.2@data$Longitud != 0 , ]

# remove data with missing altitude data
todos.4<-todos.3[todos.3@data$Altitud!=9999,] 

# remove inconsistent data samples
todos.5<-todos.4[todos.4@data$Validacion!="Inconsistente",]

# create new column for species 
todos.5$maiz<-"Zea mays mays"

# make maices object
todos.6<-todos.5


# subset races with greater than 20 samples
maices = subset(todos.6, length(todos.6$Raza_prima) > 20)
# get sample counts of races
raza.counts<-maices@data %>% group_by(Raza_prima) %>% summarise(n()) 
raza.counts

dir.create(paste0(dir_maices))
writeOGR(maices,dsn=paste0(dir_maices),layer="maices",driver="ESRI Shapefile",overwrite=TRUE)




# make presence/absence matrix
install.packages("letsR")
library(letsR)
xy<-cbind(maices@data$Longitud,maices@data$Latitud)
PA<-letsR::lets.presabip.points(xy,maices@data$Raza_prima, resol = 0.01)
plot(PA)
plot(PA$Richness_Raster)
pam<-PA$Presence_and_Absence_Matrix
View(pam)
pam[pam == 0 ] <- NA

write.csv(pam,file=paste0(dir_out,"/pam.csv"))
save.image(file=paste0(dir_out,"/clean_maices_obs.RData"))

save.image(paste0(dir_pres,"/raster_processing.RData"))


# get topographic variables
urltopo<-c("http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/")
doctopo<-htmlParse(urltopo)
#get <a> nodes.
Anodes<-getNodeSet(doctopo,"//a")

#make the full url
urls<-grep("30s",sapply(Anodes, function(Anode) xmlGetAttr(Anode,"href")),value=TRUE)

#Select the files of interest 
r<-regexpr("wc.*?zip",urls)
files<-regmatches(urls,r)
mapply(function(x,y) download.file(x,y),urls,files)
lapply(grep(".zip",files, value=TRUE),unzip)


### DOWNLOAD WORLDCLIM ALTITUDE SRTM
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



##############################################
# download and unzip present 2.0 climate data
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
                    
                    
                    # following Condo et al. 2011 "Regional climate change scenarios for MÃ©xico." AtmÃ³sfera 24(1), 125-140 (2011)
                    c("CC", # CCSM4 (Community Climate System Model, UCAR)- new version of CCSM-30 
                      "MC", # MIROC5 (Model for Interdisciplinary Research on Climate)- new version of MIROC-HI; "Although its good performance and high resolution, MIROC32-HIRES model has an inconvenience: its sensibility is 5.6 ÂºC, way higher than the 3 ÂºC marked as âbest estimateâ in IPCCÂ´s AR4 (Wigley, 2008). " (Conde et al. 2011)
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
i<-cmip5files[1]
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
  #preclabel<-strsplit(layerNames(i[[20:31]]), "[_]")[[1]][1]
  #tminlabel<-strsplit(layerNames(i[[32:43]]), "[_]")[[1]][1]
  #tmaxlabel<-strsplit(layerNames(i[[44:55]]), "[_]")[[1]][1]

  png(filename=paste0(dir_figs,"/",x,"_",biolabel,".png"))
  par(mfrow=c(4,5))
  plot(i[[1:19]])
  dev.off()
  
  #png(filename=paste0(dir_figs,"/",x,"_",preclabel,".png"))
  #par(mfrow=c(3,4))
  #plot(i[[20:31]]/10) # prcp
  #dev.off()
  
  #png(filename=paste0(dir_figs,"/",x,"_",tminlabel,".png"))
  #par(mfrow=c(3,4))
  #plot(i[[32:43]]/10) # tmin
  #dev.off()
  
  #png(filename=paste0(dir_figs,"/",x,"_",tmaxlabel,".png"))
  #par(mfrow=c(3,4))
  #plot(i[[44:55]]/10) # tmax
  #dev.off()
}


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
# BUILD MODELLING STACKS 
# write rasters to file 
pres_cropstack@title<-"1970-2000"
f50cropstack@title<-"2041-2060"
f70cropstack@title<-"2061-2080"

stackApply()
cropstacks<-c(pres_cropstack,f50cropstack,f70cropstack)
i<-cropstacks[[1]]
install.packages("envirem")
library(envirem)
f50cropstack
names(f50cropstack)
envirem::layerCreation(f50cropstack,var="all")

i<-cropstacks[[1]]
for (i in cropstacks){
  tmean<-mosaic(i[[44:55]],i[[32:43]],fun=mean)
  tmin <- i[[32:43]]
  tmax <- i[[44:55]]
  prec <- i[[20:31]]
  

  GrowingMonths<-nlayers(i[[20:31]]>100)
  names(GrowingMonths)<-"GrowingMonths"

  # calculate temperature extremes
  temp <- otherTempExtremes(tmean, tmin, tmax)
  meantempWarmest <- temp[['meanTempWarmest']]
  meantempColdest <- temp[['meanTempColdest']]
  cont<-continentality(meantempWarmest, meantempColdest)
  names(cont)<-"continentality"
  
  # growing degree days
  ggd<-growingDegDays(tmean, 10)
  names(ggd)<-"GrowingDegreeDays"
  
  alt<-raster(paste0(dir_topo,"/alt_cropped.grd"))
  names(alt)<- "elev"
  terrain<-terrain(alt, opt=c('slope','TRI'), unit='degrees', neighbors=8)
  names(terrain)
  
  p<-stack(i,GrowingMonths,cont,ggd,alt,terrain)
  names(p)<-c(layerNames(i),"GrowingMonths","cont","ggd","alt","terrain")
  writeRaster(p, paste0(dir_stacks,"/",i@title,".grd"), bylayer=FALSE, format='raster', overwrite=T)
  assign(i, p)
}

#biolabel i[[1:19]]
#prec i[[20:31]]
#tmin i[[32:43]]
#tmax i[[44:55]]



vifstep<-vifstep(pres_cropstack,th=10)

vifcor<-vifcor(pres_cropstack,th=0.75)

# remove multicollinearity of full stack
library(SpaDES)
library(virtualspecies)
coll_vars_pres<-virtualspecies::removeCollinearity(pres_cropstack,multicollinearity.cutoff = 0.45, select.variables = FALSE, sample.points = F,plot = TRUE)
save.image(paste0(dir_clim,"/raster_processing.RData"))
print(coll_vars_pres)
coll_vars_pres


# manually select raster layers to retain from climate pca, including other variables of interest
presmodstack<-raster::stack(paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_01.tif"),
                            paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_02.tif"),
                            paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_03.tif"),
                            paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_05.tif"),
                            paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_06.tif"),
                            paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_15.tif"),
                            paste0(dir_p.mosaics,"crop/crop_wc2.0_bio_30s_18.tif"),
                            paste0(dir_topo,"/alt_cropped.grd"),
                            paste0(dir_ind,"/pob-ind.grd")
                            
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




library(raster)
save(presmodstack,file=paste0(dir_stacks,"/present_modstack.RData"))
writeRaster(presmodstack, paste0(dir_stacks,"/present_modstack.grd"), bylayer=FALSE, format='raster',overwrite=TRUE)

save(f50modstack,file=paste0(dir_stacks,"/f50_modstack.RData"))
writeRaster(f50modstack, paste0(dir_stacks,"/f50_modstack.grd"), bylayer=FALSE, format='raster',overwrite=TRUE)

save(f70modstack,file=paste0(dir_stacks,"/f70_modstack.RData"))
writeRaster(f70modstack, paste0(dir_stacks,"/f70_modstack.grd"), bylayer=FALSE, format='raster',overwrite=TRUE)


save.image(file=paste0(dir_clim,"/raster_processing.RData"))











































#### OLD COPY raster stack cropping

# Raster Stack Cropping

img<-list.files(paste0(dir_pres,"/2.0"),pattern="wc", full.names=TRUE)
img

topo<-list.files(paste0(dir_pres,"/2.0"),pattern="topo", full.names=TRUE)
topo

# stack climatic variables
clim1<-stack(img)

# stack topographic variables
topo1<-stack(topo)


# read in some data for mexico to get outline
pob_ind_2010<-readOGR(dsn=paste0(dir_dat,"/indigeneity/CONABIO/Municipio/m_pob_ind_2010_poinmun10gw.shp"),layer="m_pob_ind_2010_poinmun10gw")
View(pob_ind_2010@data)
plot(pob_ind_2010)
require(rgeos)

# dissolve mexico municipio boundaries 
mex<-gUnaryUnion(pob_ind_2010)
plot(mex)

# get bounding box
box1<-as(raster::extent(mex), "SpatialPolygons")

# buffer bounding box
# box2<-gBuffer(box1, width=0.5, byid=TRUE )
# plot(box1)
# plot(box2,add=TRUE)
# since its in wgs84, corners are rounded, so get extent of box
# bbox<-as(raster::extent(box2), "SpatialPolygons")
# plot(box1)
# plot(bbox,add=TRUE)
# plot(stack1[[1]],add=TRUE)

# clip stack to buffered bounding box
clim2<-crop(clim1,box1, snap='near')
clim2<-stack(clim2)

# crop topographic variables and stack
topo2<-crop(topo1,bbox, snap='out')
topo2<-stack(topo2)

# create cropped stack
cropstack<-stack(clim2,topo2)
# plot(stack3[1])

# cropstack stack
# fullstack<-stack(clim1,topo1)


cropstack@layers
cropstack[[1:12]]

bios<-biovars(cropstack[[1:12]],cropstack[[49:60]],cropstack[[37:48]])
biostack<-stack(cropstack,bios)


save(biostack,file=paste0(dir_pres,"/2.0/crop/cropstack.R"))
save.image(paste0(dir_clim,"/raster_processing.RData"))

writeRaster(biostack, paste0(dir_pres,"/2.0/crop/cropstack.tif"), bylayer=FALSE, format='GTiff')
writeRaster(biostack, paste0(dir_pres,"/2.0/crop/crop.tif"), bylayer=TRUE, suffix = names(biostack),format='GTiff')


# need to run with more RAM
install.packages(("virtualspecies"))
require(virtualspecies)
coll_vars<-virtualspecies::removeCollinearity(biostack,multicollinearity.cutoff = 0.65,  sample.points = TRUE, plot = FALSE)
save.image(paste0(dir_clim,"/raster_processing.RData"))

print(coll_vars)

# new vars
stack3<-stack(paste0(dir_pres,"/2.0/crop/","/crop_X30as__present_tmean_mosaic_layer.05.tif"),
              paste0(dir_pres,"/2.0/crop/","wc2.0_30s_prec_05.tif"),
              paste0(dir_pres,"/2.0/crop/","/crop-_X30as__present_biovar_mosaic__bio04.tif"),
              paste0(dir_pres,"/2.0/crop/","/crop-_X30as__present_biovar_mosaic__bio05.tif"),
              paste0(dir_pres,"/2.0/crop/","/crop-_X30as__present_biovar_mosaic__bio06.tif"),
              paste0(dir_pres,"/2.0/crop/","/crop-_X30as__present_biovar_mosaic__bio08.tif"),
              paste0(dir_pres,"/2.0/crop/","/crop-_X30as__present_biovar_mosaic__bio09.tif"),
              paste0(dir_pres,"/2.0/crop/","bio03.tif"),
              paste0(dir_pres,"/2.0/crop/","bio12.tif"),
              paste0(dir_pres,"/2.0/crop/","bio17.tif"),
              paste0(dir_pres,"/2.0/crop/","bio19.tif"),
              paste0(dir_pres,"/2.0/crop/","/crop-_X30as__topo_roughness_mosaic.tif")
)

save(stack3,file=paste0(dir_out,"/cropstack-mcl-rm.R"))
writeRaster(stack3, paste0(dir_out,"/cropstack-mcl-rm.tif"), bylayer=FALSE, format='GTiff')
writeRaster(stack3, paste0(dir_out,"/cropstack-mcl-rm.tif"), bylayer=TRUE, suffix = names(stack3),format='GTiff')

###################################################################################################

# With all data extracted for each tile and variable into one folder, 
# renamed with Bulk Rename Utility to add zeros before single digital 
# month identifiers in name, so that files can be listed in sequential 
# order. The suffix of the file indicates the tile. Stack each raster 
# tile with pattern="".

tile11data=list.files(dir_p.raw, pattern="*_11.tif$", full.names=TRUE)
s.tile11<-stack(tile11data)
dim(s.tile11)

tile12data=list.files(dir_p.raw, pattern="*_12.tif$", full.names=TRUE)
s.tile12<-stack(tile12data)
dim(s.tile12)


tile13data=list.files(dir_p.raw, pattern="*_13.tif$", full.names=TRUE)
s.tile13<-stack(tile13data)
dim(s.tile13)


tile21data=list.files(dir_p.raw, pattern="*_21.tif$", full.names=TRUE)
s.tile21<-stack(tile21data)
dim(s.tile21)


tile22data=list.files(dir_p.raw, pattern="*_22.tif$", full.names=TRUE)
s.tile22<-stack(tile22data)
dim(s.tile22)


tile23data=list.files(dir_p.raw, pattern="*_23.tif$", full.names=TRUE)
s.tile23<-stack(tile23data)
s.tile23

# 48 variables for tmin, tmean, tmax and prec, corresponding to each month. Altitude is only one layer, so adds one to each raster stack.

#mosaic stacks
tile_stacks_list<-list(s.tile11,s.tile12,s.tile13,s.tile21,s.tile22,s.tile23)

tile_stacks_list$filename<-'present_tp_30acs'
tile_stacks_list$fun<-'mean'


present_tp_30acs<-do.call(merge,tile_stacks_list)
dim(present_tp_30acs)
plot(present_tp_30acs)
present_tp_30acs_stack<-stack(present_tp_30acs)
writeRaster(present_tp_30acs_stack,filename="present_",format="GTiff",bylayer=TRUE,suffix=names(present_stack_30acs))

present_p_30acs<-subset(present_stack_30acs,2:13)
plot(present_p_30acs)
present_tmin_30acs<-subset(present_stack_30acs,38:49)
present_tmax_30acs<-subset(present_stack_30acs,14:25)
#create biovars
present_bio_stack<-biovars(present_p_30acs,present_tmin_30acs,present_tmax_30acs)

plot(present_bio_stack)


#stack t and prec with biovars to create brick
present_brick_vars_30acs<-brick(present_alt_stack,present_prec_stack,present_)

```


```{r}
#alternative, longer way. different dimensions?? than other raster stack above
s.tile11_alt<-subset(s.tile11,1)
s.tile12_alt<-subset(s.tile12,1)
s.tile13_alt<-subset(s.tile13,1)
s.tile21_alt<-subset(s.tile21,1)
s.tile22_alt<-subset(s.tile22,1)
s.tile23_alt<-subset(s.tile23,1)

s.tile11_prec<-subset(s.tile11,2:13)
s.tile12_prec<-subset(s.tile12,2:13)
s.tile13_prec<-subset(s.tile13,2:13)
s.tile21_prec<-subset(s.tile21,2:13)
s.tile22_prec<-subset(s.tile22,2:13)
s.tile23_prec<-subset(s.tile23,2:13)



s.tile11_tmin<-subset(s.tile11,38:49)
dim(s.tile11_tmin)
s.tile12_tmin<-subset(s.tile12,38:49)
s.tile13_tmin<-subset(s.tile13,38:49)
s.tile21_tmin<-subset(s.tile21,38:49)
s.tile22_tmin<-subset(s.tile22,38:49)
s.tile23_tmin<-subset(s.tile23,38:49)



s.tile11_tmax<-subset(s.tile11,14:25)
dim(s.tile11_tmax)
s.tile12_tmax<-subset(s.tile12,14:25)
s.tile13_tmax<-subset(s.tile13,14:25)
s.tile21_tmax<-subset(s.tile21,14:25)
s.tile22_tmax<-subset(s.tile22,14:25)
s.tile23_tmax<-subset(s.tile23,14:25)


s.tile11_tmean<-subset(s.tile11,26:37)
s.tile12_tmean<-subset(s.tile12,26:37)
s.tile13_tmean<-subset(s.tile13,26:37)
s.tile21_tmean<-subset(s.tile21,26:37)
s.tile22_tmean<-subset(s.tile22,26:37)
s.tile23_tmean<-subset(s.tile23,26:37)


mosaic.alt<-mosaic(s.tile11_alt,s.tile12_alt,s.tile13_alt,s.tile21_alt,s.tile22_alt,s.tile23_alt,fun="mean")

mosaic.prec<-mosaic(s.tile11_prec,s.tile12_prec,s.tile13_prec,s.tile21_prec,s.tile22_prec,s.tile23_prec,fun="mean")

mosaic.tmin<-mosaic(s.tile11_tmin,s.tile12_tmin,s.tile13_tmin,s.tile21_tmin,s.tile22_tmin,s.tile23_tmin,fun="mean")

mosaic.tmax<-mosaic(s.tile11_tmax,s.tile12_tmax,s.tile13_tmax,s.tile21_tmax,s.tile22_tmax,s.tile23_tmax,fun="mean")

mosaic.tmean<-mosaic(s.tile11_tmean,s.tile12_tmean,s.tile13_tmean,s.tile21_tmean,s.tile22_tmean,s.tile23_tmean,fun="mean")

writeRaster(mosaic.alt,filename="present_altitude_mosaic",format="GTiff")
writeRaster(mosaic.prec,filename="present_prec_mosaic",format="GTiff",bylayer=TRUE,suffix=names(mosaic.prec))
writeRaster(mosaic.tmin,filename="present_tmin_mosaic",format="GTiff",bylayer=TRUE,suffix=names(mosaic.tmin))
writeRaster(mosaic.tmax,filename="present_tmax_mosaic",format="GTiff",bylayer=TRUE,suffix=names(mosaic.tmax))
writeRaster(mosaic.tmean,filename="C:/Users/Steven Gonzalez/Desktop/Geo-7300/Data/climatological/present/mosaics/present_tmean_mosaic",format="GTiff",bylayer=TRUE,suffix=names(mosaic.tmean))

#rename mosaics in bulk rename utility so suffix names are sequential; add zero before single digits.
dim(mosaic.alt)
mosaic.alt

dim(mosaic.prec)
mosaic.prec

dim(mosaic.tmin)
mosaic.tmin

dim(mosaic.tmax)
mosaic.tmax

dim(mosaic.tmean)
mosaic.tmean

#Calculate biovars
present_bio_30acs_sm<-biovars(mosaic.tmean,mosaic.prec,mosaic.tmean,mosaic.tmin,mosaic.biovars)
mosaic_bios<-stack(mosaic.tmean)
writeRaster(mosaic_bios,filename="C:/Users/Steven Gonzalez/Desktop/Geo-7300/Data/climatological/present/mosaics/present_bios_mosaic",format="GTiff",bylayer=TRUE,suffix=names(mosaic_bios))
