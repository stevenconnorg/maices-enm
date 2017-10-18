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

# 
dir_clim<-paste0(dir_dat,"/clim")
dir_pres<-paste0(dir_clim,"/present")
dir_fut<-paste0(dir_clim,"/future")

dir_p.mosaics<-paste0(dir_pres,"/2.0/")
dir_f.mosaics<-paste0(dir_fut,"/1.4/")

dir_stacks<-paste0(dir_clim,"/stacks/")

folders<-as.list(ls())

# tk work on function to create these directories
for (i in 1:length(folders))  { 
  f<-folders[[i]]
  folder<-get(f)
  dir.create(folder,recursive=TRUE) 
} 


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

#####################################
# Indigenous vars processing and rasterization using CONABIO data
# http://www.conabio.gob.mx/informacion/gis/maps/geo

setwd(dir_ind)

# define a vector of variable names wanted to download
# copy link address of .shp you want to download to view var name from CONABIO geoportal

vars<-c("poinmun10gw", # Indigenous population data by municipio geometry
        "lengmun90gw", # 1st, 2nd, 3rd, and 4th major indigenous language by municipio geometry (1990)
        "presindigw"   # Categorical indicators of indigenous population magnitude by municipio
)



for (i in vars) {
  download.file(paste0("http://www.conabio.gob.mx/informacion/gis/maps/geo/",i,".zip"),destfile = paste0(dir_ind,"/",i,".zip"),method = "wget")
  zipF<-paste0(dir_ind,"/",i,".zip") # lets you choose a file and save its file path in R (at least for windows)
  outDir<-paste0(dir_ind,"/",i) # Define the folder where the zip file should be unzipped to 
  unzip(zipF,exdir=outDir)  # unzip your file 
  
}

# for (m in replacements) {
#  files <- list.files(outDir,pattern = i, full.names = F) 
#   sapply(files,FUN=function(eachPath){
#   file.rename(from=eachPath,to=sub(pattern= paste0(i,".*"),replacement= paste0(m,".*"),eachPath))
#     })
#   }

}

ind_shps<-array()

# read in shapefiles inside dir (recursive), naming Spatial* object with shpfile name
shps <- dir(dir_ind, "*.shp",recursive = T)
shps
shps <- gsub('\\.shp$',"",shps)
for (shp in shps) assign(shp, readOGR(paste0(dir_ind,"/",shp,"/",shp,".shp"),layer=shp))


## 

#read in shapefiles
shapefiles <- dir(dir_ind, "*.shp",recursive=TRUE) # read in shapefiles into r environment


## Do pca on ethnographic variables
??pca
# read in some data for mexico to get outline to crop
library(rgeos)
library(tibble)
library(OpenStreetMap)
library(tmap)
pob_ind_2010<-readOGR(dsn=paste0(dir_dat,"/indigeneity/CONABIO/Municipio/poinmun10gw.shp"),layer="poinmun10gw")
pob_ind_2010@data$pcnt1<-(pob_ind_2010@data$P3I10 / pob_ind_2010@data$P3T10)*100
pob_ind_2010@data$pcnt2<-(pob_ind_2010@data$P3I10 / pob_ind_2010@data$P3T10)*100
pob_ind_2010@data$pcnt3<-(pob_ind_2010@data$MON10 / pob_ind_2010@data$P3T10)*100
pob_ind_2010@data$pcnt4<-(pob_ind_2010@data$BIL10 / pob_ind_2010@data$P3T10)*100
pob_ind_2010@data$pcnt5<-(pob_ind_2010@data$NEI10 / pob_ind_2010@data$P3T10)*100
pob_ind_2010@data$pcnt6<-(pob_ind_2010@data$NLI10 / pob_ind_2010@data$P3T10)*100
pob_ind_2010@data$pcnt7<-(pob_ind_2010@data$NELI10 / pob_ind_2010@data$P3T10)*100
pca<-prcomp(~pcnt1+pcnt2+pcnt3+pcnt4+pcnt5+pcnt6+pcnt7,data= pob_ind_2010@data,center=TRUE,scale=TRUE)
print(pca)
plot(pca,type="l")
screeplot(pca)
summary(pca)

# plot pca values from dataframe
library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)
g <- ggbiplot(pca, obs.scale = 1, var.scale = 1, 
              # groups = ir.species, 
              ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

# join pca to original shape
pob_ind_2010@data<-rownames_to_column(pob_ind_2010@data,var="id")
scores<-as.data.frame(pca$x)
scores<-rownames_to_column(scores,var="id")

# merge scores to spolydf
pob_ind_2010@data<-merge(x=pob_ind_2010@data,y=scores,by="id")
colnames(pob_ind_2010@data)

# plot pc1 over geometry
qtm(shp = pob_ind_2010, fill = c(paste0("PC1",rep(1:7))), fill.palette = "-Blues",ncol=2) 
# crs(pob_ind_2010)
# ind_wgs = spTransform(pob_ind_2010, CRS("+init=epsg:4326"))
ind_wgs<-pob_ind_2010
osm_tiles = tmaptools::read_osm(bbox(ind_wgs)) #

# plot pc1 over geometry with osm tiles
tm_shape(osm_tiles) + tm_raster() + tm_shape(ind_wgs) +
  tm_fill("PC1", fill.title = "Ethnolinguistic Component", scale = 0.8, alpha = 0.5) +
  tm_layout(legend.position = c(0.89,0.02))


# rasterize pca results for pc 1 and 2
ethnlang1 <- rasterize(pob_ind_2010, r, field = pob_ind_2010@data$PC1, fun = "mean", 
                       update = FALSE, updateValue = "NA")
ethnlang2 <- rasterize(pob_ind_2010, r, field = pob_ind_2010@data$PC2, fun = "mean", 
                       update = FALSE, updateValue = "NA")


# plot rasters
plot(ethnlang1)
plot(ethnlang2)

extent(cropstack)
extent(ind)
extent(ind) <- extent(cropstack)

writeRaster(ethnlang1, paste0(dir_dat,"/ethn/rasts/ethn_lang-pc1.tif"), bylayer=FALSE, format='GTiff',overwrite=TRUE)
plot(ind)

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

# http://www.worldclim.org/tiles.php
altrasters<-Sys.glob(file.path("C:/Users/Steven Gonzalez/Desktop/Geo-7300/Data/climatological/present/raw/tiffs/alt*.tif"))
alt.list<-list()
for(i in 1:length(altrasters)) {alt.list[i]<-raster(altrasters[i])}

alt.list$fun <- mean
alt.mosaic <- do.call(mosaic,alt.list)

writeRaster(alt.mosaic,filename="C:/Users/Steven Gonzalez/Desktop/Geo-7300/Data/climatological/present/mosaics/altitude_mosaic.tif",format="GTiff")/
  plot(alt.mosaic,main="Altitude (m) (30 arc seconds x 30 arc seconds")

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
writeRaster(prescropstack, paste0(dir_stacks,"/present_cropstack.grd"), bylayer=FALSE, format='raster',overwrite=F)

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
                  model=c("HE","IP","HD","HG"), # models
                  year=c(50,70), ##50 or 70     # period 2050 or 2070
                  res="30s") %>%                # resolution
  
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
for(i in zfs[17:32]){ 
  exdir= gsub(".zip","",i)
  unzip(i,exdir=exdir)  # unzip file
  # unlink(i)  # remove zip file
  patt<-substr(i,nchar(i)-12+1,nchar(i)-4)
  gtifs<-list.files(exdir,pattern=patt,full.names=T)# [c(2, 6:13, 3:5)]##reorder
  tempstack<-stack(grds) ##rasterbrick
  ctempstack<-crop(tempstack,bbox) ##filtered
  writeRaster(ctempstack,bylayer=TRUE,filename=paste0(dir_f.mosaics,"/crop/",patt,".tif"),overwrite=TRUE)
  # unlink(gtifs)
  print(paste0("Finished with file ",patt," (",which(zfs==i)," out of ",length(zfs),")"))
  do.call(file.remove,list(list.files(pattern="temp*"))) 
  for (i in gtifs[i]){
    writeRaster(i,format="raster",filename=paste0(dir_f.mosaics,"/",patt,".grd")
  }
  grds<-list.files(paste0(dir_f.mosaics,pattern="grd",full.names=T)
                   unlink(gtifs)
}


# ensemble GCMs by monthly means
# get list of all tifs in crop directory, mean ensemble by variable, and write layer to stack
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
  for (t in topofiles){
    file.copy(paste0(dir_dat,"/topo/",t),paste0(dir_f.mosaics,"/crop/ensemble/",period,"/",basename(t)))
  }
}


#  make present crop raster stack
f50_croprasts<-list.files(paste0(dir_f.mosaics,"crop/ensemble/50/"),pattern="\\.tif$", full.names=TRUE)
f70_croprasts<-list.files(paste0(dir_f.mosaics,"/crop/ensemble/70/"),pattern="\\.tif$", full.names=TRUE)

f50cropstack<-stack(f50_croprasts)
f70cropstack<-stack(f70_croprasts)

# write to file
save(f50cropstack,file=paste0(dir_stacks,"/f50cropstack.RData"))
save(f70cropstack,file=paste0(dir_stacks,"/f70cropstack.RData"))
hdr(f50cropstack, format = "ENVI") # to preserve layer names in other programs with .grd/.gri file types -- uncompressed, so they take a while
hdr(f70cropstack, format = "ENVI")
writeRaster(f50cropstack, paste0(dir_stacks,"/f50cropstack.grd"), bylayer=FALSE, format='raster', overwrite=F)
writeRaster(f70cropstack, paste0(dir_stacks,"/f70cropstack.grd"), bylayer=FALSE, format='raster', overwrite=F)
save.image(paste0(dir_clim,"/raster_processing.RData"))


#####################################
# remove multicollinearity of full stack
# need to run with more RAM
grds<-list.files(dir_stacks,pattern=".grd")

f50cropstack<-stack(paste0(dir_stacks,"/",grds[1]))
f70cropstack<-stack(paste0(dir_stacks,"/",grds[2]))
pres_cropstack<-stack(paste0(dir_stacks,"/",grds[3]))

library(SpaDES)
layerNames(present_cropstack)

install.packages(("virtualspecies"))
library(virtualspecies)
coll_vars_pres<-virtualspecies::removeCollinearity(pres_cropstack,multicollinearity.cutoff = 0.7, select.variables = FALSE, sample.points = TRUE, nb.points = (pres_cropstack@ncols/2),plot = TRUE)
save.image(paste0(dir_clim,"/raster_processing.RData"))
print(coll_vars_pres)
coll_vars_pres
# manually select raster layers to retain from climate pca, including other variables of interest
presmodstack<-stack(paste0(dir_p.mosaics,"/crop/crop_bio5.tif"),
                    paste0(dir_p.mosaics,"/crop/crop_bio6.tif"),
                    paste0(dir_p.mosaics,"/crop/crop_bio18.tif"),
                    paste0(dir_p.mosaics,"/crop/crop_bio3.tif"),
                    paste0(dir_p.mosaics,"/crop/crop_bio17.tif"),
                    paste0(dir_p.mosaics,"/crop/crop_bio4.tif"),
                    paste0(dir_p.mosaics,"/crop/crop_bio15.tif"),
                    paste0(dir_p.mosaics,"/crop/crop_bio8.tif"),
                    paste0(dir_p.mosaics,"/crop/crop_bio9.tif"),
                    paste0(dir_p.mosaics,"/crop/crop_bio2.tif"),
                    paste0(dir_dat,"/topo/crop/crop_topo_roughness_mosaic.tif"),
                    paste0(dir_dat,"/topo/crop/crop_topo_altitude_mosaic.tif"),
                    paste0(dir_dat,"/ethn/rasts/ethn_lang-pc1.tif")
                    
)

coll_vars_pres_mod<-virtualspecies::removeCollinearity(presmodstack,multicollinearity.cutoff = 0.7, select.variables = FALSE, sample.points = TRUE, nb.points = (presmodstack@ncols/2),plot = TRUE)

presmodstack

f50modstack<-stack(paste0(dir_f.mosaics,"crop/ensemble/50/85bi50_5_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/50/85bi50_6_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/50/85bi50_18_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/50/85bi50_3_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/50/85bi50_17_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/50/85bi50_4_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/50/85bi50_15_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/50/85bi50_15_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/50/85bi50_8_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/50/85bi50_9_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/50/85bi50_2_ensemble.tif"),
                   paste0(dir_dat,"/topo/crop/crop_topo_roughness_mosaic.tif"),
                   paste0(dir_dat,"/topo/crop/crop_topo_altitude_mosaic.tif"),
                   paste0(dir_dat,"/ethn/rasts/ethn_lang-pc1.tif")
                   
                   
)

f70modstack<-stack(paste0(dir_f.mosaics,"crop/ensemble/70/85bi70_5_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/70/85bi70_6_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/70/85bi70_18_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/70/85bi70_3_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/70/85bi70_17_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/70/85bi70_4_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/70/85bi70_15_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/70/85bi70_8_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/70/85bi70_9_ensemble.tif"),
                   paste0(dir_f.mosaics,"crop/ensemble/70/85bi70_2_ensemble.tif"),
                   paste0(dir_dat,"/topo/crop/crop_topo_roughness_mosaic.tif"),
                   paste0(dir_dat,"/topo/crop/crop_topo_altitude_mosaic.tif"),
                   paste0(dir_dat,"/ethn/rasts/ethn_lang-pc1.tif")
)


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
