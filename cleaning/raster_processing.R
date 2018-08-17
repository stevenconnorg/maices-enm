library(raster)
library(XML)
library(rgdal)
# define directories

root <- "E:/thesis"
#root <- "~/thesis"

setwd(root)

dir_dat <- paste0(root, "/01_data")
dir_R <- paste0(root, "/02_R")
dir_out <- paste0(root, "/03_output")
dir_figs <- paste0(root, "/04_figs")
dir_lit <- paste0(root, "/05_lit")
dir_comp <- paste0(root, "/06_comp")
dir_presentations <- paste0(root, "/07_pres")

dir_maices <- paste0(dir_dat, "/maices")
dir_ind <- paste0(dir_dat, "/ind")
dir_topo <- paste0(dir_dat, "/topo")
#
dir_clim <- paste0(dir_dat, "/clim")
dir_pres <- paste0(dir_clim, "/present")
dir_fut <- paste0(dir_clim, "/future")

dir_p.mosaics <- paste0(dir_pres, "/2.0/")
dir_f.mosaics <- paste0(dir_fut, "/1.4/")

dir_land <- paste0(dir_dat, "/land-cover")

dir_stacks <- paste0(dir_dat, "/stacks/")
dir_envirem <- paste0(dir_dat, "/envirem/")

folders <- as.list(ls())
i <- folders[[3]]
# tk work on function to create these directories
for (i in 1:length(folders))  {
  f <- folders[[i]]
  folder <- get(f)
  dir.create(folder, recursive = TRUE)
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
r <- raster(nrow = 2304, ncol = 3770)
extent(r) <- c(-117.625, -86.20833, 14.025, 33.225)
r


##############################################


##############################################
# get soil data from FAO harmonilzed world soil db v1.2
##############################################

urlsoil <-
  c(
    "http://www.fao.org/soils-portal/soil-survey/soil-maps-and-databases/harmonized-world-soil-database-v12/en/"
  )
docsoil <- htmlParse(urlsoil)
#get <a> nodes.
Anodes <- getNodeSet(docsoil, "//a")
grep("*.asc*", Anodes, value = T)
#make the full url
urls <-
  grep("asc", sapply(Anodes, function(Anode)
    xmlGetAttr(Anode, "href")), value = TRUE)
urls <- paste0("http://www.fao.org/", urls)
filename <- basename(urls)
library(tools)
ext <- file_ext(filename)
asc <- paste0(".", ext[1])
file <- gsub(asc, "", filename)
landfiles <- paste0(dir_land, "/", filename)
mapply(function(x, y)
  download.file(x, y), urls, files)

r <- raster(nrow = 2304, ncol = 3770)
extent(r) <- c(-117.625, -86.20833, 14.025, 33.225)


bbox <- bbox(r)
for (i in landfiles) {
  ##140:173
  name <- gsub(".asc", "", basename(i))
  ras <- raster(i)
  ri <- resample(ras, r)
  ri <- round(ri, 0)
  croplay <- crop(ri, bbox) ##filtered
  writeRaster(
    croplay,
    filename = paste0(dir_land, "/crop/crop_", name, ".grd"),
    overwrite = TRUE
  )
  do.call(file.remove, list(list.files(pattern = "temp*")))
}

# rename and reorder stack layers

landfiles <-
  list.files(paste0(dir_land, "/crop"),
             full.names = T,
             pattern = ".grd")
landstack <- stack(landfiles)
names(landstack)

names(landstack) <-
  c(
    "TotalCult",
    "IrrCult",
    "Rain-fedCult",
    "Forested",
    "Grass/Woodland",
    "Barren",
    "SQ1NutAvail",
    "SQ2NutRetn",
    "SQ3RootingCond",
    "SQ4OxyAvailRoots",
    "SQ5ExcessSalts",
    "SQ6Toxicity",
    "SQ7Workability",
    "Urban",
    "Water"
  )
#plot(landstack)
landcoverstack <- landstack[[c(1:6, 14:15)]] #extract land cover layers
writeRaster(
  landcoverstack,
  filename = paste0(dir_stacks, "/FAOlandstack.grd"),
  bylayer = FALSE,
  format = "raster",
  overwrite = TRUE
)

soilqualstack <- landstack[[7:13]] # extract soil quality variables
library(raster)
for (i in 1:nlayers(soilqualstack)) {
  soilqualstack[[i]] <-
    raster::as.factor(soilqualstack[[i]]) # assign soil quality variables as factors
  
}
writeRaster(
  soilqualstack,
  filename = paste0(dir_stacks, "/FAOsoilstack.grd"),
  bylayer = FALSE,
  format = "raster",
  overwrite = TRUE
)

##############################################
# Download WorldClim altitude Rasters by tile
##############################################
tiles <- (c("00","01","02","03","04","05","06","07","08","09","010","011",
          "10","11", "12", "13", 
                                  "14","15","16","17","18","19","110","111",
          "20",
                "21", "22", "23",
                                  "24","25","26","27","28","29","210","211",
          "30","31","32","33","34","35","36","37","38","39","310","311",
          "40","41","42","43","44","45","46","47","48","49","410","411"
))

for (i in tiles) {
  download.file(
    paste0(
      "http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/alt_",
      i,
      "_tif.zip"
    ),
    destfile = paste0(dir_topo, "/alt_", i, ".zip"),
    method = "wget"
  )
  zipF <-
    paste0(dir_topo, "/alt_", i, ".zip") # lets you choose a file and save its file path in R (at least for windows)
  outDir <-
    paste0(dir_topo, "/alt_", i) # Define the folder where the zip file should be unzipped to
  unzip(zipF, exdir = outDir)  # unzip your file
}

# altrasters<-list.files(dir_topo,pattern=".bil",full.names = T,recursive=T)

altrasters <-
  lapply(list.files(
    dir_topo,
    pattern = ".tif",
    full.names = T,
    recursive = T
  ),
  raster)

altrasters$fun <- mean
alt.mosaic <- do.call(mosaic, altrasters)
crs(alt.mosaic) <-
  "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
names(alt.mosaic) <- "altitude"

# alt.crop <- crop(alt.mosaic, r)
# writeRaster(
#   alt.crop,
#   filename = paste0(dir_topo, "/alt_cropped.grd"),
#   format = "raster",
#   overwrite = T
# )

writeRaster(
  alt.mosaic,
  filename = paste0(dir_topo, "/alt.grd"),
  format = "raster",
  overwrite = T
)

#alt.crop <- raster(paste0(dir_topo, "/alt_cropped.grd"))
alt.crop <- raster(paste0(dir_topo, "/alt_cropped.grd"))

names(alt.crop) <- "elev"
terrain <-
  terrain(
    alt.crop,
    opt = c('aspect', 'slope', 'Roughness'),
    unit = 'degrees',
    neighbors = 8
  )

## transform aspect to linear feature with cosine (Northness)
asp.radians <- terrain[['aspect']] * pi / 180
plot(asp.radians)
cos.asp.rad <- cos(asp.radians)
plot(cos.asp.rad)
names(cos.asp.rad) <- "cos.asp.rad"
writeRaster(
  cos.asp.rad,
  paste0(dir_stacks, "cos.asp.rad.grd"),
  format = 'raster',
  overwrite = T
)

#Sappington et al., (2007) vector ruggedness measure library(spatialEco)
vrm <- spatialEco::vrm(alt.crop)
writeRaster(
  vrm,
  paste0(dir_stacks, "vrm.grd"),
  bylayer = FALSE,
  format = 'raster',
  overwrite = T
)

names(vrm) <- "vrm"

topostack <-
  stack(terrain[['roughness']], terrain[['slope']], cos.asp.rad, vrm)

writeRaster(
  topostack,
  paste0(dir_stacks, "topostack-full.grd"),
  bylayer = FALSE,
  format = 'raster',
  overwrite = T
)



##############################################
# download and unzip present Wordclim 2.0 climate data
##############################################
setwd(dir_p.mosaics)
url <- c("http://worldclim.org/version2")
doc <- htmlParse(url)
#get <a> nodes.
Anodes <- getNodeSet(doc, "//a")

#make the full url
urls <-
  grep("30s", sapply(Anodes, function(Anode)
    xmlGetAttr(Anode, "href")), value = TRUE)

#Select the files of interest
r <- regexpr("wc.*?zip", urls)
files <- regmatches(urls, r)
mapply(function(x, y)
  download.file(x, y), urls, files)
lapply(grep(".zip", files, value = TRUE), unzip)

# alternatively with a cluster
# get inside present clim data mosaic directory
# setwd(dir_p.mosaics)

# library(parallel)
# cl <- parallel::makeCluster(detectCores())
# clusterMap(cl, download.file, url = urls, destfile = paste0(dir_p.mosaics,"/",files))

# stopCluster(cl)

#####
# Present clim data Raster Stack Cropping
dir.create(paste0(dir_p.mosaics, "/crop"))

# read in present data, stack and write to RData file
pres_rasts <-
  list.files(paste0(dir_p.mosaics),
             pattern = "\\.tif$",
             full.names = TRUE)

presfullstack <- stack(pres_rasts)

# write to file
save(presfullstack, file = paste0(dir_stacks, "/present_fullstack.RData"))
writeRaster(
  presfullstack,
  paste0(dir_stacks, "/present_fullstack.grd"),
  bylayer = FALSE,
  format = 'GTiff',
  overwrite = FALSE
)







bbox <- extent(r)
# i<-pres_rasts[1]
pres_rasts <-
  list.files(paste0(dir_p.mosaics),
             pattern = "srad",
             full.names = TRUE)

for (i in pres_rasts) {
  ##140:173
  name <- gsub(".tif", "", basename(i))
  ras <- raster(i)
  croplay <- crop(ras, bbox) ##filtered
  writeRaster(
    croplay,
    filename = paste0(dir_p.mosaics, "/crop/crop_", name, ".tif"),
    overwrite = TRUE
  )
  do.call(file.remove, list(list.files(pattern = "temp*")))
}

save.image(paste0(dir_R, "/raster_processing.RData"))

#  make present crop raster stack
pres_croprasts <-
  list.files(paste0(dir_p.mosaics, "/crop"),
             pattern = "\\.tif$",
             full.names = TRUE)

pres_cropstack <- stack(pres_croprasts)


# write to file
save(pres_cropstack, file = paste0(dir_stacks, "/full_present_cropstack.RData"))
writeRaster(
  prescropstack,
  paste0(dir_stacks, "/present_cropstack.grd"),
  bylayer = FALSE,
  format = 'raster',
  overwrite = T
)

bio <-
  pres_cropstack[[c(grep(
    layerNames(pres_cropstack),
    pattern = "bio",
    value = T
  ))]]

tmin <-
  pres_cropstack[[c(grep(
    layerNames(pres_cropstack),
    pattern = "tmin",
    value = T
  ))]]
tmax <-
  pres_cropstack[[(grep(
    layerNames(pres_cropstack),
    pattern = "tmax",
    value = T
  ))]]
prec <-
  pres_cropstack[[(grep(
    layerNames(pres_cropstack),
    pattern = "prec",
    value = T
  ))]]

#tmin<-tmin*10
#tmax<-tmax*10

pres_cropstack<-stack(bio,prec,tmin,tmax)


#divideLayers<-function(layers){
#for (i in layers){
#    cropstack[[i]] <-cropstack[[i]]*10
#   
#}

save(pres_cropstack, file = paste0(dir_stacks, "/present_cropstack.RData"))
writeRaster(
  pres_cropstack,
  paste0(dir_stacks, "/present_cropstack.grd"),
  bylayer = FALSE,
  format = 'raster',
  overwrite = T
)


#### DOWNLOAD FUTURE DATA FOR RCP 85 ####
# from https://github.com/ClimateActionUCI/datasets/blob/master/get_CMIP5.R

library(dplyr)
library(raster)
AR5 <- 'AR5 temperature increase projections'
ar5.df <-
  data.frame(
    'Scenario' = c('RCP2.6', 'RCP4.5', 'RCP6.0', 'RCP8.5'),
    '2046 - 2065' = c(
      '1.0 (0.4 to 1.6)',
      '1.4 (0.9 to 2.0)',
      '1.3 (0.8 to 1.8)',
      '2.0 (1.4 to 2.6)'
    ),
    '2081 - 2100' = c(
      '1.0 (0.3 to 1.7)',
      '1.8 (1.1 to 2.6)',
      '1.2 (1.4 to 3.1)',
      '3.7 (2.6 to 4.8)'
    )
  )
mods <-
  expand.grid(
    var = c("tn", "tx", "pr", "bi"),
    #tn, tx, pr, or bi, no bi?
    rcp = c(85),
    ##26, 45, 60, 85   # rcp
    model =
      
      
      # following Condo et al. 2011 "Regional climate change scenarios for MÃÂ©xico." AtmÃÂ³sfera 24(1), 125-140 (2011)
      c("CC", # CCSM4 (Community Climate System Model, UCAR)- new version of CCSM-30
        "MC", # MIROC5 (Model for Interdisciplinary Research on Climate)- new version of MIROC-HI; "Although its good performance and high resolution, MIROC32-HIRES model has an inconvenience: its sensibility is 5.6 ÃÂºC, way higher than the 3 ÃÂºC marked as Ã¢ÂÂbest estimateÃ¢ÂÂ in IPCCÃÂ´s AR4 (Wigley, 2008). " (Conde et al. 2011)
        "MP", # MPI-ESM-LR (Max-Plank Institute) - per 5th National Communication of Mexico for the United Nations Framework Convention on Climate Change http://atlasclimatico.unam.mx/atlas/Docs/f_escenarios.html
        "HE", # HADGEM2-ES (Met Office Hadley)per 5th  removed because already downloaded
        "GF") # GFDL-CM3 (Geophysical Fluid Dynamics Laboratory )
    ,
    
    year = c(50, 70),
    ##50 or 70     # period 2050 or 2070
    res = "30s"
  )                # resolution

mutate(
  filename = paste0(tolower(model), rcp, var, year),
  url = paste0(
    "http://biogeo.ucdavis.edu/data/climate/cmip5/",
    res,
    "/",
    filename,
    ".zip"
  )
)

# set path of where to download future data zips to
outpath <- dir_f.mosaics

# set download.file timeout
options(timeout = 10000)
# make function to download files
dwnldfxn <- function(aurl, filename) {
  try(raster:::.download(aurl, filename))
}

# extract urls from mods df
urls <- mods$url

# get zipfiles
zipfile <- paste0(outpath, substr(urls, nchar(urls) - 12 + 1, nchar(urls)))


# download each zipfile for each url
mapply(dwnldfxn, aurl = urls, filename = zipfile)

###########################################################
## unzip files and crop
###########################################################

# for cropping, read in some data for mexico to get outline
bbox <- extent(r)

# subset
path = paste0(dir_f.mosaics)
zfs <- list.files(path, pattern = "zip", full.names = T)
str(zfs)

dir.create(paste0(dir_f.mosaics, "/crop/ensemble/50"), recursive = TRUE)
dir.create(paste0(dir_f.mosaics, "/crop/ensemble/70"), recursive = TRUE)

# unzip zip directories for each climatology into respective directories, crop, and write layers to 'crop' directory
for (i in zfs) {
  i<-zfs[1]
  ##140:173
  exdir = gsub(".zip", "", i)
  unzip(i, exdir = exdir)  # unzip file
  apatt <- substr(i, nchar(i) - 12 + 1, nchar(i) - 4)
  gtifs <-
    list.files(exdir, pattern = ".tif", full.names = T)# [c(2, 6:13, 3:5)]##reorder
  for (t in gtifs) {
    gtifras <- raster(t)
    cropped <- crop(gtifras, bbox)
    writeRaster(
      cropped,
      format = "raster",
      filename = paste0(dir_f.mosaics, "crop/", gsub(".tif", "", basename(t)), ".grd"),
      overwrite = T
    )
    # unlink(t)
    do.call(file.remove, list(list.files(pattern = "temp*")))
  }
  do.call(file.remove, list(list.files(pattern = "temp*")))
}




# ensemble GCMs by monthly means
# get list of all tifs in crop directory, mean ensemble by variable, and write layer to stack
cmip5files <-
  list.files(paste0(dir_f.mosaics, "/crop"),
             pattern = "*.grd",
             recursive = FALSE)

#for (i in patterns) {

ensembleRas<- function(patterns){
  pattern<-patterns[1]
  i<-pattern
  
  
    # clims <- gsub(".grd", "", basename(i))
    varfiles1 <-
      list.files(paste0(dir_f.mosaics, "/crop"),
                 pattern = paste0(pattern,"\\.grd"),
                 recursive = F)
    paste0("Using files :")
    varfiles1 
    varfiles2 <- paste0(dir_f.mosaics, "crop/", varfiles1)
    name <- gsub(".grd", "", basename(varfiles2))
    name <- paste0(substr(name, 3, nchar(name)))
    name<-name[1]
    
    
    stack <- stack(varfiles2)
    paste0("Raster stack consists of :")
    stack
    paste0("Averaging raster stack for :")
    name
    
    stack2 <-
      stackApply(stack, indices = rep(1:1, length(names(stack))), fun = "mean")
    names(stack2) <- paste0(name, "_ensemble")

    period <- substr(i, 3,4)
    paste0("Period :")
    period
    writeRaster(
      stack2,
      filename = paste0(
        dir_f.mosaics,
        "crop/ensemble/",
        period,
        "/",
        name,
        "_ensemble.grd"
      ),
      overwrite = TRUE
    )
}
library(snowfall)
library(parallel)
sfInit( parallel=TRUE, cpus=detectCores()-1)

# # Export packages to snowfall
print("Exporting packages to cluster")
sfLibrary('raster', character.only = TRUE)

## do monthly mean variables first (1:19)
patterns<-c("pr50","pr70","tn50","tn70","tx70","tx50")
months<-c("1","2","3","4","5","6","7","8","9","10","11","12")
vars.df<-expand.grid(patterns,months)
patterns<-paste0(vars.df$Var1,vars.df$Var2)
patterns
sfExport("dir_f.mosaics","patterns","ensembleRas")


sfLapply(patterns,ensembleRas)

# then do for bio vars  (1:19)
patterns<-c("bi50","bi70")
months<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")
vars.df<-expand.grid(patterns,months)
patterns<-paste0(vars.df$Var1,vars.df$Var2)

sfExport("dir_f.mosaics","patterns","ensembleRas")

sfLapply(patterns,ensembleRas)

sfStop()

#  make  crop raster stack
f50_croprasts <-
  list.files(
    paste0(dir_f.mosaics, "crop/ensemble/50/"),
    pattern = "\\.grd$",
    full.names = TRUE
  )
f70_croprasts <-
  list.files(
    paste0(dir_f.mosaics, "/crop/ensemble/70/"),
    pattern = "\\.grd$",
    full.names = TRUE
  )

f50cropstack <- stack(f50_croprasts)
f70cropstack <- stack(f70_croprasts)

# write to file
# 
# save(f50cropstack, file = paste0(dir_stacks, "/f50cropstack.RData"))
# save(f70cropstack, file = paste0(dir_stacks, "/f70cropstack.RData"))
# hdr(f50cropstack, format = "ENVI") # to preserve layer names in other programs with .grd/.gri file types -- uncompressed, so they take a while
# hdr(f70cropstack, format = "ENVI")
# writeRaster(
#   f50cropstack,
#   paste0(dir_stacks, "/f50cropstack.grd"),
#   bylayer = FALSE,
#   format = 'raster',
#   overwrite = T
# )
# writeRaster(
#   f70cropstack,
#   paste0(dir_stacks, "/f70cropstack.grd"),
#   bylayer = FALSE,
#   format = 'raster',
#   overwrite = T
# )


### READ IN ALL RASTER STACKS, OVERWRITE WITH MATCHING LAYER INDICES
library(raster)
library(quickPlot)
library(usdm)
# read in raster stacks

 grds <- list.files(path = dir_stacks,
                     pattern = ".grd",
                     full.names = T)
  
# pres_cropstack <- stack(grds[19])
# f50cropstack <- stack(grds[3])
# f70cropstack <- stack(grds[6])

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

# if you need to subset a raster stack to match another
    # preslayervecs <- c(73:91, 1:12, 37:48, 25:36)
    # subset(pres_cropstack)
    # layerNames(pres_cropstack[[preslayervecs]])
    # nlayers(pres_cropstack[[preslayervecs]])
    # pres_cropstack <- pres_cropstack[[preslayervecs]]
    # layerNames(pres_cropstack)

# make sure raster stack names are the same, formatting first
library(quickPlot)
names(pres_cropstack) <-
  gsub("crop_wc2.0", "", layerNames(pres_cropstack)) # remove prefix
names(pres_cropstack) <-
  gsub("30s_", "", layerNames(pres_cropstack))        # remove suffix
names(pres_cropstack) <-
  gsub("X_", "", layerNames(pres_cropstack))        # remove suffix
layerNames(pres_cropstack)

names(f50cropstack)
names(f70cropstack)
#names(f50cropstack[[c(1,12:19,2:11,20,24:31,21:23,32,36:43,33:35,44,48:55,45:47)]])
f50cropstack<-f50cropstack[[c(1,12:19,2:11,20,24:31,21:23,32,36:43,33:35,44,48:55,45:47)]]
names(f50cropstack)

f70cropstack<-f70cropstack[[c(1,12:19,2:11,20,24:31,21:23,32,36:43,33:35,44,48:55,45:47)]]
names(f70cropstack)

layerNames(pres_cropstack)

names(f50cropstack) <- layerNames(pres_cropstack)  # apply presmodstack layer names to future stacks
names(f70cropstack) <- layerNames(pres_cropstack)

# save raster stacks
f50brick <- brick(f50cropstack)
f70brick <- brick(f70cropstack)
presbrick <- brick(pres_cropstack)

hdr(f50brick, format = "ENVI") # to preserve layer names in other programs with .grd/.gri file types -- uncompressed, so they take a while
hdr(f70brick, format = "ENVI")
hdr(presbrick, format = "ENVI")

save(pres_cropstack,file=paste0(dir_stacks, "/present_cropstack.RData"))
save(f50cropstack,file=paste0(dir_stacks, "/f50cropstack.RData"))
save(f70cropstack,file=paste0(dir_stacks, "/f70cropstack.RData"))

# writeRaster(
#   presbrick,
#   paste0(dir_stacks, "/present_cropstack.grd"),
#   bylayer = FALSE,
#   format = 'raster',
#   overwrite = T
# )
# writeRaster(
#   f50brick,
#   paste0(dir_stacks, "/f50cropstack.grd"),
#   bylayer = FALSE,
#   format = 'raster',
#   overwrite = T
# )
# writeRaster(
#   f70brick,
#   paste0(dir_stacks, "/f70cropstack.grd"),
#   bylayer = FALSE,
#   format = 'raster',
#   overwrite = T
# )
save.image(paste0(dir_clim, "/raster_processing.RData"))


### MAKE SOME PLOTS
library(raster)
# read in raster stacks
# grds <- list.files(path = dir_stacks,
#                    pattern = ".grd",
#                    full.names = T)
# grds
# pres_cropstack <- stack(grds[19])
# f50cropstack <- stack(grds[11])
# f70cropstack <- stack(grds[16])

# write rasters to file
pres_cropstack@title <- "1970-2000"
f50cropstack@title <- "2041-2060"
f70cropstack@title <- "2061-2080"
names(f70cropstack)
fcropstacks <- list(f50cropstack, f70cropstack)
library(rasterVis)
library(quickPlot)


writeRaster(
  f50cropstack,
  paste0(dir_stacks, "/f50cropstack.grd"),
  bylayer = FALSE,
  format = 'raster',
  overwrite = T
)
writeRaster(
  f70cropstack,
  paste0(dir_stacks, "/f70cropstack.grd"),
  bylayer = FALSE,
  format = 'raster',
  overwrite = T
)

save(f50cropstack,file=paste0(dir_stacks, "/f50cropstack.RData"))
save(f70cropstack,file=paste0(dir_stacks, "/f70cropstack.RData"))
library(raster)
f50cropstack<-stack(paste0(dir_stacks, "/f50cropstack.grd"))
f70cropstack<-stack(paste0(dir_stacks, "/f70cropstack.grd"))

# or load RData
load(file=paste0(dir_stacks, "/present_cropstack.RData"))
load(file=paste0(dir_stacks, "/f50cropstack.RData"))
load(file=paste0(dir_stacks, "/f70cropstack.RData"))

# layerRanges <- list(c(1:19), c(20:31), c(32:43), c(44:55))
# 
# for (i in cropstacks[i]) {
#   for (l in layerRanges[[l]]) {
#     as.integer(l)
#     layers <- as.character(l)
#     x <- i@title
#     label <- strsplit(layerNames(i[[l]]), "[_]")[[1]][1]
#     png(filename = paste0(dir_figs, "/", x, "_", label, ".png"))
#     plot(i[[l]]) # prcp
#     dev.off()
#   }
# }



#fcropstacks

                # for (cropstack in fcropstacks) {
                #   preclabel <- strsplit(layerNames(cropstack[[20:31]]), "[_]")[[1]][1]
                #   tminlabel <- strsplit(layerNames(cropstack[[32:43]]), "[_]")[[1]][1]
                #   tmaxlabel <- strsplit(layerNames(cropstack[[44:55]]), "[_]")[[1]][1]
                #   
                #   # get max temp diff
                #   tmaxdiff <- ((f70cropstack[[44:55]]) - ((pres_cropstack[[44:55]])))
                #   png(
                #     filename = paste0(
                #       dir_figs,
                #       "/",
                #       tmaxlabel,
                #       " (C) 1970 - 2000 avg. to ",
                #       cropstack@title,
                #       " avg. Difference, RCP 8.5.png"
                #     )
                #   )
                #   rasterVis::levelplot(
                #     tmaxdiff,
                #     main = paste0(
                #       tmaxlabel,
                #       " (C) 1970 - 2000 avg. to ",
                #       cropstack@title,
                #       " avg. Difference, RCP 8.5"
                #     )
                #   )
                #   dev.off()
                #   
                #   # get max temp diff
                #   tmindiff <- (cropstack[[32:43]] ) - (pres_cropstack[[32:43]])
                #   png(
                #     filename = paste0(
                #       dir_figs,
                #       "/",
                #       tminlabel,
                #       " (C) 1970 - 2000 avg. to ",
                #       cropstack@title,
                #       " avg. Difference, RCP 8.5.png"
                #     )
                #   )
                #   rasterVis::levelplot(
                #     tmindiff,
                #     main = paste0(
                #       tminlabel,
                #       " (C) 1970 - 2000 avg. to ",
                #       cropstack@title,
                #       " avg. Difference, RCP 8.5"
                #     )
                #   )
                #   dev.off()
                #   
                #   # get prcp diff
                #   pcpdiff <- (cropstack[[20:31]]) - (pres_cropstack[[20:31]])
                #   png(
                #     filename = paste0(
                #       dir_figs,
                #       "/",
                #       pclabel,
                #       " (mm) 1970 - 2000 avg. to ",
                #       cropstack@title,
                #       " avg. Difference, RCP 8.5.png"
                #     )
                #   )
                #   rasterVis::levelplot(
                #     pcpdiff,
                #     main = paste0(
                #       pclabel,
                #       " 1970 - 2000 avg. to ",
                #       cropstack@title,
                #       " avg. Difference, RCP 8.5"
                #     )
                #   )
                #   dev.off()
                #   
                # }

#####################################
# Generate envirem data
# write rasters to file

library(raster)
# grds <-
#   list.files(
#     path = paste0(dir_pres, "/2.0/crop"),
#     pattern = ".tif",
#     full.names = T
#   )
# 
# presbio_cropstack <- stack(grep(grds, pattern = "bio", value = T))
# presprec_cropstack <- stack(grep(grds, pattern = "prec", value = T))
# prestmin_cropstack <- stack(grep(grds, pattern = "tmin", value = T))
# prestmax_cropstack <- stack(grep(grds, pattern = "tmax", value = T))


# 
# pres_cropstack <-
#   stack(presbio_cropstack,
#         presprec_cropstack,
#         prestmin_cropstack,
#         prestmax_cropstack)
# 
# 

# library(quickPlot)
# names(pres_cropstack) <-
#   gsub("crop_wc2.0", "", layerNames(pres_cropstack)) # remove prefix
# names(pres_cropstack) <-
#   gsub("30s_", "", layerNames(pres_cropstack))        # remove suffix
# names(pres_cropstack) <-
#   gsub("X_", "", layerNames(pres_cropstack))        # remove suffix
# layerNames(pres_cropstack)


load(file=paste0(dir_stacks, "/f50cropstack.RData"))
load(file=paste0(dir_stacks, "/f70cropstack.RData"))

pres_rasts <-
  list.files(paste0(dir_p.mosaics),
             pattern = "\\.tif$",
             full.names = TRUE)

presfullstack <- stack(pres_rasts)


# layerNames(f50cropstack)
# layerNames(f70cropstack)

pres_cropstack@title <- "1970-2000"
f50cropstack@title <- "2041-2060"
f70cropstack@title <- "2061-2080"
presfullstack@title<-"1970-2000"

library(envirem)
library(quickPlot)
library(stringr)
pres_cropstack


## take full present stack and rewrite layers to tiff to fit envirem data




for (i in list(presfullstack,pres_cropstack)) {
  #presfullstack <-stack(list.files(paste0(dir_stacks,"1970-2000/global"),full.names = T))
  f50stack <-stack(list.files(paste0(dir_stacks,"2041-2060"),full.names = T))
  presfullstack[[1:19]]
  i<-f50stack
  bio <-
    i[[(grep(
      layerNames(i),
      pattern = "bio",
      value = T
    ))]]
  tmin <-
    i[[(grep(
      layerNames(i),
      pattern = "tmin",
      value = T
    ))]]
  tmax <-
    i[[(grep(
      layerNames(i),
      pattern = "tmax",
      value = T
    ))]]
  prec <-
    i[[(grep(
      layerNames(i),
      pattern = "prec",
      value = T
    ))]]
  
  srad <-
    i[[(grep(
      layerNames(i),
      pattern = "srad",
      value = T
    ))]]
  
  
  tmean <-
    i[[(grep(
      layerNames(i),
      pattern = "tavg",
      value = T
    ))]]
  
  ### the present datalayers are WorldClim 2.0, with Temp vars in C * 10
  ## multiply to match future WorldClim 1.4 data
  
  bio[[c(1,2,4:11)]]<-bio[[c(1,2,4:11)]]*10
  tmax<-tmax*10
  tmean<-tmean*10
  tmin<-tmin*10
  
  p <- stack(bio, prec, srad, tmax, tmean, tmin)
  
  
  p@title <- i@title
  names(p) <- c(
    paste0("bio_", rep(01:19)),
    paste0("prec_", rep(01:12)),
    paste0("et_solrad_", rep(01:12)),
    paste0("tmax_", rep(01:12)),
    paste0("tmean_", rep(01:12)),
    paste0("tmin_", rep(01:12))
  )
  
  
  dir <- paste0(dir_stacks, i@title,"/global")
  dir.create(dir, showWarnings = FALSE)
  writeRaster(
    p,
    paste0(dir, "/", names(p), ".tif"),
    bylayer = TRUE,
    format = 'GTiff',
    overwrite = T
  )
  
  
}


for (i in list(presfullstack,pres_cropstack)) {
  
  dir <- paste0(dir_stacks, i@title,"/global")
  
  rasterExt = '.tif'
  files <- list.files(path = dir, pattern = paste0(rasterExt, '$'))


  # https://stat.ethz.ch/pipermail/r-sig-geo/2011-May/011793.html
  rasterOptions(maxmemory=3e+10, chunksize=1e+07)
  
  var <- c('minTempWarmest',
           'monthCountByTemp10',
           'PETColdestQuarter',
           'PETDriestQuarter',
           'PETWarmestQuarter')
  
  generateRasters(
    var = var,
    maindir = dir,
    resName = NULL,
    timeName = i@title,
    outputDir = paste0(dir_dat, "/envirem/global"),
    rasterExt = ".tif",
    nTiles = 1,
    overwriteResults = TRUE,
    outputFormat = "GTiff",
    tempDir = "~/temp",
    gdalinfoPath = NULL,
    gdal_translatePath = NULL
  )

}


# do the same for the  cropped stacks
for (i in list(f50cropstack,f70cropstack)) {
  bio <-
    i[[(grep(
      layerNames(i),
      pattern = "bio",
      value = T
    ))]]
  tmin <-
    i[[(grep(
      layerNames(i),
      pattern = "tmin",
      value = T
    ))]]
  tmax <-
    i[[(grep(
      layerNames(i),
      pattern = "tmax",
      value = T
    ))]]
  prec <-
    i[[(grep(
      layerNames(i),
      pattern = "prec",
      value = T
    ))]]
  
  tmean<-tmin
    for (l in 1:nlayers(tmin)){
    tmean[[l]] <- mosaic(tmin[[l]], tmax[[l]], fun = mean)
  }

  names(tmean) <-
    paste0("tmean_", str_pad(rep(01:12), 2, pad = "0"))
  tmean
  
  sradFiles <-
    list.files(paste0(dir_pres, "/2.0/crop"),
               pattern = "srad",
               full.names = T)
  srad <- stack(sradFiles)
  names(srad) <- paste0("srad_", str_pad(rep(01:12), 2, pad = "0"))
  
  
  p <- stack(bio, prec, srad, tmax, tmean, tmin)
 

  p@title <- i@title
  names(p) <- c(
    paste0("bio_", rep(01:19)),
    paste0("prec_", rep(01:12)),
    paste0("et_solrad_", rep(01:12)),
    paste0("tmax_", rep(01:12)),
    paste0("tmean_", rep(01:12)),
    paste0("tmin_", rep(01:12))
  )
  
  
  dir <- paste0(dir_stacks, i@title)
  dir.create(dir, showWarnings = FALSE)
  
  writeRaster(
    p,
    paste0(dir, "/", names(p), ".tif"),
    bylayer = TRUE,
    format = 'GTiff',
    overwrite = T
  )
  
    
  rasterExt = '.tif'
    	files <- list.files(path = dir, pattern = paste0(rasterExt, '$'))
  	
  # bioclimFiles <- grep('bio_\\d\\d?', files, value = TRUE)
  # 	bioclimFiles <- gsub('(bio_\\d\\d?)(\\.\\w+$)', '\\1', bioclimFiles)
  # 	if (all(paste0('bio_', 1:19) %in% bioclimFiles) & length(bioclimFiles) == 19) {
  # 		bioclimCheck <- TRUE
  # 	} else {
  # 		bioclimCheck <- FALSE
  # 	}
    	
    	
  # https://stat.ethz.ch/pipermail/r-sig-geo/2011-May/011793.html
  rasterOptions(maxmemory=3e+07, chunksize=1e+06)
  
  var <- c('minTempWarmest',
    'monthCountByTemp10',
    'PETColdestQuarter',
    'PETDriestQuarter' ,
    'PETWarmestQuarter')
  
    generateRasters(
      var = var,
      maindir = dir,
      resName = NULL,
      timeName = ,i@title,
      outputDir = paste0(dir_dat, "/envirem/"),
      rasterExt = ".tif",
      nTiles = 1,
      overwriteResults = TRUE,
      outputFormat = "GTiff",
      tempDir = "~/temp",
      gdalinfoPath = NULL,
      gdal_translatePath = NULL
    )
}



library(raster)

# read in new cropstacks with all calculated variables
f70cropstack@title

for (cropstack in list(f50cropstack, f70cropstack)) {
 
   print(cropstack@title)
  
 biostack <-
    raster::stack(list.files(
      path = paste0(dir_stacks, cropstack@title, "/"),
      pattern = "bio",
      full.names = T
    ))
 
  enviremstack <-
    raster::stack(list.files(
     path = paste0(dir_envirem, cropstack@title, "/"),
      pattern = ".tif",
      full.names = T
    ))
  
  ### CLEAN UP LAYER NAMES AFTER ENVIREM ###
  # get biostack names back to equal lengths
  names(biostack)[nchar(names(biostack)) == 5] <-
    sub("(.{4})(.*)", "\\10\\2", layerNames(biostack)[nchar(layerNames(biostack)) ==
                                                        5])
  names(biostack)
  # and reorder
  biostack <- biostack[[c(1, 12:19, 2:11)]]
  
  names(biostack)
  names(enviremstack)
  
  full_stack <- stack(biostack, enviremstack)
  
  names(full_stack)
  # change envirem time id to match cropstack title with hyphen
  time <- gsub("-", ".", cropstack@title)
  # remove id's added by envirem
  names(full_stack) <- gsub(paste0("X", time, "__")
                            , ""
                            , names(full_stack))
  
  # set 19 bioclim vars with real name
  names(full_stack)     <- c(
    "BIO1AnnMeanTemp",
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
    "EVMminTempWarmest",
    "EVMmonthCountByT10",
    "EVMPETColdestQ",
    "EVMPETDriestQ",
    "EVMPETWarmestQ"
  )
  
  writeRaster(
    full_stack,
    paste0(dir_stacks, "/", cropstack@title, "_full-biostack", ".grd"),
    bylayer = FALSE,
    format = 'raster',
    overwrite = T
  )
  
}


full_bio_names     <- c(
    "BIO1AnnMeanTemp",
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

bionames     <- c(
    "BIO1AnnMeanTemp",
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
    "EVMminTempWarmest",
    "EVMmonthCountByT10",
    "EVMPETColdestQ",
    "EVMPETDriestQ",
    "EVMPETWarmestQ"
  )

###########################################################
###       MODELLING STACK VARIABLE SELECTION            ###
###########################################################

### GET full present (cropped) modelling stack
# bio - WorldClim and ENVIREM layers
# FAO land stack
# topographic variables

pres_biostack <- stack(paste0(dir_stacks, "/1970-2000biostack.grd"))
landstack <- stack(paste0(dir_stacks, "/FAOlandstack.grd"))
topostack <- stack(paste0(dir_stacks, "/topostack.grd"))

names(pres_biostack)

fullpresstack <- stack(pres_biostack, landstack, topostack)

# create correlation matrix of bioclimatic variables
library(corrplot)
fullpresstack[is.na(fullpresstack)] <- 0

# OR try
funNA <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}
fullpresstack <- calc(fullpresstack, funNA)

# OR maybe
fullpresstack <- reclassify(fullpresstack, cbind(NA, 0))



mat <- as.matrix(fullpresstack)
colnames(mat) <- names(fullpresstack)

#View(mat)
cormat <- cor(mat)
# visualize correlation matrix

png(
  filename = paste0(dir_figs, "/biovars_cormat.png"),
  width = 1000,
  height = 1000
)
par(mfrow = c(1, 1))
corrplot(cormat)
dev.off()



write.csv(mat, file = paste0(dir_out, "/fullpresstack_matrix.csv"))
#mat<-read.csv(file=paste0(dir_out,"/fullpresstack_matrix.csv"))
write.csv(cormat, file = paste0(dir_out, "/fullpresstack_corr_matrix.csv"))

###########################################################
###       MAXENT VARIABLE SELECTION (UNDER CONSTRUCTION)          
###########################################################

maxentvarselect<-function(){

##### DO MAXENT VARIABLE SELECTION #####
library(MaxentVariableSelection)

# read in all observations across varieties

library(biomod2)
library(rgdal)

maices<-readOGR(dsn=paste0(dir_maices,"/todos-maices-cleaned.shp"),layer="todos-maices-cleaned")
maiz.enm<- enmtools.species()

equalarea<-CRS("+proj=aea +lat_1=14.5 +lat_2=32.5 +lat_0=24 +lon_0=-105 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
latlong<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

maices_EA<-spTransform(maices,equalarea)
xy_EA<-coordinates(maices_EA)
head(maices@data$maiz)
# use template wgs raster for environmental predictors to create presence/absence matrix

r<-raster(nrow=2304,ncol=3770)
extent(r)<-c(-117.625,-86.20833,14.025,33.225)
r
crs(r) <- latlong
# project raster to equal area
rSp<-projectRaster(r,crs=equalarea)


# get extent of object
rSp_extent<-extent(rSp) 

# spatial filtering of points at 5 times the resolution < 5km^2
PA_EA<-letsR::lets.presab.points(xy_EA,maices@data$maiz,xmn= rSp_extent@xmin, xmx= rSp_extent@xmax   , ymn=rSp_extent@ymin , ymx= rSp_extent@ymax, resol = res(rSp)*5,crs =equalarea,show.matrix=TRUE,remove.cells=TRUE,remove.sp=TRUE)
XYobs <-PA_EA[,1:2]
xymaices<-SpatialPoints(coordinates(XYobs))
crs(xymaices)<-equalarea
latlongmaices<-spTransform(xymaices,latlong)
latlongobs<-coordinates(latlongmaices)


#library(raster)
#Grids <- fullpresstack 
# Extracting the variables for all occurrencelocations
#VariablesAtOccurrencelocations <- raster::extract(projectRaster(Grids,rSp),XYobs)

# Combining the extracted values with the longitude and latitude values
#Outfile <- as.data.frame(cbind("Zeamays", latlongobs,
#VariablesAtOccurrencelocations))
#colnames(Outfile) <- c("species","longitude","latitude",
#colnames(VariablesAtOccurrencelocations))
#Outfile <-Outfile [complete.cases(Outfile),]

#write.csv(Outfile, file = paste0(root,"/02_R/maices-enm/varselect/VariablesAtOccurrencelocations.csv"), append = FALSE,sep = ",", eol = "\n", na = "NA", dec = ".",col.names = TRUE,row.names=FALSE)

#library(dismo)
#mask <- raster(fullpresstack[[1]])

# select 500 random points
# set seed to assure that the examples will always
# have the same random sample.
#bg <- randomPoints(mask, 100000 )

#mexico<-readOGR(dsn=paste0(dir_dat,"/mexico-outline.shp"),layer="mexico-outline")

# make spatial points from bg locations
#spbg<-SpatialPoints(bg)

# apply crs
#crs(spbg)<-crs(mexico)

# subset points inside mexico outline
#spbg<-spbg[mexico] 

# extract variables at point locations
#VariablesAtBGlocations <- raster::extract(Grids,coordinates(spbg))

# bind data with species name and lat/long coordinates
#bgOutfile <- as.data.frame(cbind("Background", coordinates(spbg),
#VariablesAtBGlocations))

# rename columns to fix maxent
#colnames(bgOutfile ) <- c("species","longitude","latitude",
#colnames(VariablesAtOccurrencelocations))

# remove any points that have NA from a layer
#bgOutfile <-bgOutfile [complete.cases(bgOutfile),]
#View(bgOutfile)
# write to file
#write.csv(bgOutfile , file = paste0(root,"/02_R/maices-enm/varselect/VariablesAtBGlocations.csv"), col.names = TRUE,row.names=FALSE)


#install.packages("MaxentVariableSelection")
#library(MaxentVariableSelection)
#VariableSelection(maxent= paste0(root,"/02_R/maices-enm/maxent.jar"),
# outdir =  paste0(root,"/02_R/maices-enm/varselect/out"),
# gridfolder = gridfolder,
#occurrencelocations = paste0(root,"/02_R/maices-enm/varselect/VariablesAtOccurrencelocations.csv"),
# backgroundlocations=paste0(root,"/02_R/maices-enm/varselect/VariablesAtBGlocations.csv"),
# additionalargs="noproduct",
#contributionthreshold = 0.5,
# correlationthreshold = 0.75,
# betamultiplier = c(1, 2.5, 5, 10)
#)


## select relevant variables from AICc selection

}


###########################################################
###       VARIANCE INFLATION FACTOR CHECK               ###
###########################################################

# do final VIF check

library(usdm)
# get variance inflation factors of remaining variables
presvif<-usdm::vif(fullpresstack)
save(presvif,file=paste0(dir_out,"/vif.RData"))
load(paste0(dir_out, "/vif.RData"))
paste0(presvifstep@results$Variables)


###########################################################
###       FINAL VARIABLES TO BE USED IN MODELING        ###
###########################################################

vars <- c(
  'BIO8MeanTWettestQ',
  'BIO13PrecWettestMonth',
  'BIO15PrecSeas.COV',
  'BIO18PrecWarmestQ',
  'BIO19PrecColdestQ',
  'EVMminTempWarmest',
  'EVMmonthCountByT10',
  'EVMPETColdestQ',
  'EVMPETDriestQ' ,
  'EVMPETWarmestQ',
  'IrrCult' ,
  'Rain.fedCult'     ,
  'Grass.Woodland',
  'Barren',
  'Urban',
  'Water',
  'cos.asp.rad',
  'slope',
  'vrm'
)



###########################################################
###               GET MODELLING STACKS                  ###
###########################################################



### FOR PRESENT CROPPED STACK ###

# manually select raster layers to retain from climate pca, including other variables of interest
#presmodstack<-fullpresstack[[paste0(presvifstep@results$Variables)]]
presmodstack <- fullpresstack[[vars]]
names(presmodstack)

presmodstack <- stack(presmodstack,
                           paste0(dir_ind, "/pob-ind.grd"))




### FOR PRESENT UNCROPPED STACK ###

### get uncropped modeling stack
full.biostack <-
  raster::stack(list.files(
    path = paste0(dir_stacks,"/1970-2000/global"),
    pattern = "bio",
    full.names = T
  ))

enviremstack@title="2061-2080"
time <- gsub("-", ".", enviremstack @title)

full.biostack<-full.biostack[[c(1,12:19,2:11)]]

enviremstack <-
  raster::stack(list.files(path = paste0(dir_envirem, "/global/1970-2000/"),
    pattern = ".tif",
    full.names = T
  ))

pres_full-biostack <- stack(biostack,enviremstack)



landstack <- stack(paste0(dir_stacks, "/FAOlandstack.grd"))


full.landstack<-stack(paste0(dir_dat,"/land-cover/CULT_2000.asc"),
                      paste0(dir_dat,"/land-cover/CULTIR_2000.asc"),
                      paste0(dir_dat,"/land-cover/CULTRF_2000.asc"),
                      paste0(dir_dat,"/land-cover/FOR_2000.asc"),
                      paste0(dir_dat,"/land-cover/GRS_2000.asc"),
                      paste0(dir_dat,"/land-cover/NVG_2000.asc"),
                      paste0(dir_dat,"/land-cover/URB_2000.asc"),
                      paste0(dir_dat,"/land-cover/WAT_2000.asc")
)

full.landstack<-resample(full.landstack,full.biostack,method="bilinear")

names(full.landstack)<-names(landstack)
crs(full.landstack)<-crs(landstack)
plot(full.landstack)

writeRaster(full.landstack,filename = paste0(dir_stacks, "/FAOlandstack-full.grd"))
landstack <- stack(paste0(dir_stacks, "/FAOlandstack.grd"))
topostack <- stack(paste0(dir_stacks, "/topostack.grd"))


fullpresstack <- stack(pres_full-biostack, landstack, topostack)
pres_globmodstack <- fullf50stack[[vars]]

pres_globmodstack <-
  raster::stack(pres_globmodstack, paste0(dir_ind, "/pob-ind.grd"))




### FOR 2041 -2060  STACK ###


 biostack <-
    raster::stack(list.files(
      path = paste0(dir_stacks, "2041-2060", "/"),
      pattern = "bio",
      full.names = T
    ))

biostack<-biostack[[c(1,12:19,2:11)]]
enviremstack <-
    raster::stack(list.files(
      path = paste0(dir_envirem, "2041-2060", "/"),
      pattern = ".tif",
      full.names = T
    ))
enviremstack@title="2041-2060"
 


 time <- gsub("-", ".", enviremstack @title)
  # remove id's added by envirem
  names(enviremstack ) <- gsub(paste0("X", time, "__")
                            , ""
                            , names(enviremstack ))
  
f50biostack <- stack(biostack ,enviremstack )



names(f50biostack)<-bionames
f50landstack <- stack(paste0(dir_stacks, "/FAOlandstack.grd"))
f50topostack <- stack(paste0(dir_stacks, "/topostack.grd"))


fullf50stack <- stack(f50biostack, f50landstack, f50topostack)
names(fullf50stack)
f50modstack <- fullf50stack[[vars]]

f50modstack <-
  raster::stack(f50modstack, paste0(dir_ind, "/pob-ind.grd"))





### FOR 2061 -2080  STACK ###

 biostack <-
    raster::stack(list.files(
      path = paste0(dir_stacks, "2061-2080", "/"),
      pattern = "bio",
      full.names = T
    ))

biostack<-biostack[[c(1,12:19,2:11)]]

enviremstack <-
    raster::stack(list.files(
      path = paste0(dir_envirem, "2061-2080", "/"),
      pattern = ".tif",
      full.names = T
    ))

enviremstack@title="2061-2080"
  time <- gsub("-", ".", enviremstack @title)
  # remove id's added by envirem
  names(enviremstack ) <- gsub(paste0("X", time, "__")
                            , ""
                            , names(enviremstack ))
  
f70biostack <- stack(biostack ,enviremstack )


names(f70biostack)<-bionames

f70landstack <- stack(paste0(dir_stacks, "/FAOlandstack.grd"))
f70topostack <- stack(paste0(dir_stacks, "/topostack.grd"))


fullf70stack <- stack(f70biostack, f70landstack, f70topostack)
f70modstack <- fullf70stack[[vars]]

f70modstack <-
  raster::stack(f70modstack, paste0(dir_ind, "/pob-ind.grd"))
names(presmodstack)



library(raster)
save(presmodstack_full, file = paste0(dir_stacks, "/present_full-modstack.RData"))
writeRaster(
  presmodstack_full,
  paste0(dir_stacks, "/present_global-modstack.grd"),
  bylayer = FALSE,
  format = 'raster',
  overwrite = TRUE
)

save(presmodstack, file = paste0(dir_stacks, "/present_modstack.RData"))
writeRaster(
  presmodstack,
  paste0(dir_stacks, "/present_modstack.grd"),
  bylayer = FALSE,
  format = 'raster',
  overwrite = TRUE
)

save(f50modstack, file = paste0(dir_stacks, "/f50_modstack.RData"))
writeRaster(
  f50modstack,
  paste0(dir_stacks, "/f50_modstack.grd"),
  bylayer = FALSE,
  format = 'raster',
  overwrite = TRUE
)

save(f70modstack, file = paste0(dir_stacks, "/f70_modstack.RData"))
writeRaster(
  f70modstack,
  paste0(dir_stacks, "/f70_modstack.grd"),
  bylayer = FALSE,
  format = 'raster',
  overwrite = TRUE
)


### Project Rasters to Equal Area Projection ###
# http://spatialreference.org/ref/sr-org/38/
equalarea <-
  CRS(
    "+proj=aea +lat_1=14.5 +lat_2=32.5 +lat_0=24 +lon_0=-105 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  )

# presmodstack<-stack(paste0(dir_stacks,"present_modstack.grd"))
# f50modstack<-stack(paste0(dir_stacks,"f50_modstack.grd"))
# f70modstack<-stack(paste0(dir_stacks,"f70_modstack.grd"))

presmodstackEA <-
  projectRaster(presmodstack, crs = equalarea, method = "bilinear")
f50modstackEA <-
  projectRaster(f50modstack, crs = equalarea, method = "bilinear")
f70modstackEA <-
  projectRaster(f70modstack, crs = equalarea, method = "bilinear")


presmodstack[[1]]
writeRaster(
  presmodstackEA,
  paste0(dir_stacks, "/present_modstack_EA.grd"),
  bylayer = FALSE,
  format = 'raster',
  overwrite = TRUE
)
writeRaster(
  f50modstackEA,
  paste0(dir_stacks, "/f50_modstack_EA.grd"),
  bylayer = FALSE,
  format = 'raster',
  overwrite = TRUE
)
writeRaster(
  f70modstackEA,
  paste0(dir_stacks, "/f70_modstack_EA.grd"),
  bylayer = FALSE,
  format = 'raster',
  overwrite = TRUE
)


### Resample Raster Stacks to square pixels
#presmodstack_EA<-stack(paste0(dir_stacks,"present_modstack_EA.grd"))


# create template raster to resample to

# make sample size as original stack
r<-raster(nrow=nrow(presmodstackEA),ncol=ncol(presmodstackEA))

# same extent
extent(r)<-extent(presmodstackEA)

# but change to square resolution (1km^2)
res(r)<-1000

# apply sample crs
crs(r)<-crs(presmodstackEA)

## use bilinear if your predictors are continuous
presmodstackEA1km<-resample(presmodstackEA,r,method="bilinear")
f50modstackEA1km<-resample(f50modstackEA,r,method="bilinear")
f70modstackEA1km<-resample(f70modstackEA,r,method="bilinear")

## function to define the intersect of rasters
## remove cells that dont have data in all layers
intersect_mask <- function(x){
values_x <- raster::getValues(x)
inter_x <- values_x %*% rep(1,nlayers(x))
mask <- raster::setValues(subset(x,1),values = (inter_x>0))
return(mask)
}


## keep only all cells that are defined for all layers
print("on current stack:")

presmodstackEA1km<- stack(mask(presmodstackEA1km, intersect_mask(presmodstackEA1km)))
print("on future 50:")

f50modstackEA1km<- stack(raster::mask(f50modstackEA1km, intersect_mask(f50modstackEA1km)))
print("on future 70 stack:")

f70modstackEA1km<- stack(raster::mask(f70modstackEA1km, intersect_mask(f70modstackEA1km)))

writeRaster( presmodstackEA1km,paste0(dir_stacks, "/present_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
writeRaster(f50modstackEA1km,paste0(dir_stacks, "/f50_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
writeRaster(f70modstackEA1km,paste0(dir_stacks, "/f70_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
library(raster)

#print("Reassigning stack names:")

#names(presmodstack)<-names
#names(f50modstack)<-names
#names(f70modstack)<-names
#writeRaster( presmodstack_proj,paste0(dir_stacks, "/present_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
#writeRaster(f50modstack_proj,paste0(dir_stacks, "/f50_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
#writeRaster(f70modstack_proj,paste0(dir_stacks, "/f70_modstack_EA1km.grd"),bylayer = FALSE,format = 'raster',overwrite = TRUE)
library(raster)

presmodstack<-stack(paste0(dir_stacks,"present_modstack_EA1km.grd"))
f50modstack<-stack(paste0(dir_stacks,"f50_modstack_EA1km.grd"))
f70modstack<-stack(paste0(dir_stacks,"f70_modstack_EA1km.grd"))
print("Getting stack names:")

#names<-names(presmodstack)
names
print("Running intersect_mask function:")


#save.image(file=paste0(dir_clim,"/raster_processing.RData"))
