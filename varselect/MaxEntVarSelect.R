
install.packages("viridisLite")
install_github("danlwarren/ENMTools")
library(ENMTools)

root<-"E:/thesis"
setwd(root)
getwd()
# new directories for biomod


print("Establishing directories")
dir_dat<-paste0(root,"/01_data")
dir_R<-paste0(root,"/02_R")
dir_out<-paste0(root,"/03_output")
dir_figs<-paste0(root,"/04_figs")
dir_lit<-paste0(root,"/05_lit")
dir_comp<-paste0(root,"/06_comp")
dir_presentations<-paste0(root,"/07_pres")

dir_maices<-paste0(dir_dat,"/maices")
dir_ind<-paste0(dir_dat,"/ind")


dir_bm<-paste0(dir_R,"/00_biomod")
dir_topo<-paste0(dir_dat,"/topo")

# 
dir_clim<-paste0(dir_dat,"/clim")
dir_pres<-paste0(dir_clim,"/present")
dir_fut<-paste0(dir_clim,"/future")

dir_p.mosaics<-paste0(dir_pres,"/2.0/")
dir_f.mosaics<-paste0(dir_fut,"/1.4/")

dir_stacks<-paste0(dir_dat,"/stacks/")

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

PA_EA<-letsR::lets.presab.points(xy_EA,maices@data$maiz,xmn= rSp_extent@xmin, xmx= rSp_extent@xmax   , ymn=rSp_extent@ymin , ymx= rSp_extent@ymax, resol = res(rSp)*5,crs =equalarea,show.matrix=TRUE,remove.cells=TRUE,remove.sp=TRUE)
XYobs <-PA_EA[,1:2]
xymaices<-SpatialPoints(coordinates(XYobs))
crs(xymaices)<-equalarea
latlongmaices<-spTransform(xymaices,latlong)
latlongobs<-coordinates(latlongmaices)

pres_biostack <- stack(paste0(dir_stacks, "/1970-2000_biostack.grd"))
landstack <- stack(paste0(dir_stacks, "/FAOlandstack.grd"))
topostack <- stack(paste0(dir_stacks, "/topostack.grd"))
library(raster)

Grids <- stack(pres_biostack, landstack, topostack)

# project to equal area
GridsEA<-projectRaster(Grids,crs=equalarea)

# check crs and plot
GridsEA
plot(GridsEA)
# Extracting the variables for all occurrencelocations

VariablesAtOccurrencelocations <- raster::extract(GridsEA,XYobs)

# Combining the extracted values with the longitude and latitude values
Outfile <- as.data.frame(cbind("Zeamays", latlongobs,
VariablesAtOccurrencelocations))
colnames(Outfile) <- c("species","longitude","latitude",
colnames(VariablesAtOccurrencelocations))
Outfile <-Outfile [complete.cases(Outfile),]
View(Outfile)
write.csv(Outfile, file = paste0(root,"/02_R/maices-enm/varselect/VariablesAtOccurrencelocations.csv"), append = FALSE,sep = ",", eol = "\n", na = "NA", dec = ".",col.names = TRUE,row.names=FALSE)

library(dismo)
mask <- raster(fullpresstack[[1]])

# select 500 random points
# set seed to assure that the examples will always
# have the same random sample.
bg <- randomPoints(mask, 100000 )

mexico<-readOGR(dsn=paste0(dir_dat,"/mexico-outline.shp"),layer="mexico-outline")

# make spatial points from bg locations
spbg<-SpatialPoints(bg)

# apply crs
crs(spbg)<-crs(mexico)

# subset points inside mexico outline
spbg<-spbg[mexico] 

# extract variables at point locations
VariablesAtBGlocations <- raster::extract(Grids,coordinates(spbg))

# bind data with species name and lat/long coordinates
bgOutfile <- as.data.frame(cbind("Background", coordinates(spbg),
VariablesAtBGlocations))

# rename columns to fix maxent
colnames(bgOutfile ) <- c("species","longitude","latitude",
colnames(VariablesAtOccurrencelocations))

# remove any points that have NA from a layer
bgOutfile <-bgOutfile [complete.cases(bgOutfile),]
View(bgOutfile)
# write to file
write.csv(bgOutfile , file = paste0(root,"/02_R/maices-enm/varselect/VariablesAtBGlocations.csv"), col.names = TRUE,row.names=FALSE)


install.packages("MaxentVariableSelection")
library(MaxentVariableSelection)
VariableSelection(maxent= paste0(root,"/02_R/maices-enm/maxent.jar"),
 outdir =  paste0(root,"/02_R/maices-enm/varselect/out"),
 gridfolder = gridfolder,
occurrencelocations = paste0(root,"/02_R/maices-enm/varselect/VariablesAtOccurrencelocations.csv"),
 backgroundlocations=paste0(root,"/02_R/maices-enm/varselect/VariablesAtBGlocations.csv"),
 additionalargs="noproduct",
contributionthreshold = 40,
 correlationthreshold = 0.8,
 betamultiplier = c(2.5)
)

library(usdm)
names(fullpresstack)
vars<-names(fullpresstack)[c(5,6,13,21,25,37,38,42,46,47)]
vifstep<-vifstep(fullpresstack[[vars]],th=10)
vifstep@variables

ahli$species.name <- "ahli"
ahli$presence.points <- LonLatData 
ahli$range <- background.raster.buffer(ahli$presence.points, 50000, mask = env)
ahli$background.points <- background.points.buffer(points = ahli$presence.points,
                                                   radius = 20000, n = 1000, mask = env[[1]])