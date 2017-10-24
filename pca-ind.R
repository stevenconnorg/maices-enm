#####################################
# Indigenous vars processing and rasterization using CONABIO data
# http://www.conabio.gob.mx/informacion/gis/maps/geo

setwd(dir_ind)
library(rgdal)




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

# for (m in replacements) {  ## not working nor desrired yet
#  files <- list.files(outDir,pattern = i, full.names = F) 
#   sapply(files,FUN=function(eachPath){
#   file.rename(from=eachPath,to=sub(pattern= paste0(i,".*"),replacement= paste0(m,".*"),eachPath))
#     })
#   }

}

ind_shps<-array()

# read in shapefiles inside dir (recursive), naming Spatial* object with shpfile name
shps <- dir(dir_ind, "*.shp",recursive = T)
names <- gsub('\\.shp$',"",basename(shps))
for (name in names) assign(name, readOGR(paste0(dir_ind,"/",name,"/",name,".shp"),layer=name)) # read in sp objects with filename


## Do pca on ethnographic variables
??pca
# read in some data for mexico to get outline to crop
library(rgeos)
library(tibble)
library(OpenStreetMap)
library(tmap)

# open html metadata file with rstudio viewer
viewer <- getOption("viewer")
viewer(paste0(dir_ind,"/poinmun10gw/poinmun10gw.html"))

poinmun10gw@data$pcP3I10<-(poinmun10gw@data$P3I10 / poinmun10gw@data$P3T10)*100
# poinmun10gw@data$pcMON10<-(poinmun10gw@data$MON10 / poinmun10gw@data$P3T10)*100
# poinmun10gw@data$pcBIL10<-(poinmun10gw@data$BIL10 / poinmun10gw@data$P3T10)*100
# poinmun10gw@data$pcNEI10<-(poinmun10gw@data$NEI10 / poinmun10gw@data$P3T10)*100
# poinmun10gw@data$pcNLI10<-(poinmun10gw@data$NLI10 / poinmun10gw@data$P3T10)*100
# poinmun10gw@data$pcNELI10<-(poinmun10gw@data$NELI10 / poinmun10gw@data$P3T10)*100
# pca<-prcomp(~pcP3I10+pcP3I10+pcMON10+pcBIL10+pcNEI10+pcNLI10+pcNELI10,data= poinmun10gw@data,center=TRUE,scale=TRUE)
# print(pca)
# plot(pca,type="l")
# screeplot(pca)
# summary(pca)

# plot pca values from dataframe
library(devtools)
install_github("ggbiplot", "vqv")

# join pca to original shape
# poinmun10gw@data<-rownames_to_column(poinmun10gw@data,var="id")
# scores<-as.data.frame(pca$x)
# scores<-rownames_to_column(scores,var="id")

# merge scores to spolydf
# poinmun10gw@data<-merge(x=poinmun10gw@data,y=scores,by="id")
# colnames(poinmun10gw@data)

# plot pc1 over geometry
qtm(shp = poinmun10gw, fill = "pcP3I10", fill.palette = "-Blues",ncol=2) 

# crs(poinmun10gw)
# ind_wgs = spTransform(poinmun10gw, CRS("+init=epsg:4326"))
# ind_wgs<-poinmun10gw
osm_tiles = tmaptools::read_osm(bbox(poinmun10gw)) #

# plot pc1 over geometry with osm tiles
tm_shape(osm_tiles) + tm_raster() + tm_shape(ind_wgs) +
  tm_fill("PC1", fill.title = "Ethnolinguistic Component", scale = 0.8, alpha = 0.5) +
  tm_layout(legend.position = c(0.89,0.02))


# rasterize pca results for pc 1 and 2
ethnlang1 <- rasterize(poinmun10gw, r, field = poinmun10gw@data$PC1, fun = "mean", 
                       update = FALSE, updateValue = "NA")

popindpct <- rasterize(poinmun10gw, r, field = poinmun10gw@data$pcP3I10, fun = "mean", 
                       update = FALSE, updateValue = "NA")


# plot rasters
plot(ethnlang1)
plot(popindpct)

extent(cropstack)
extent(ind)
extent(ind) <- extent(cropstack)
names(popindpct)<-"ind_pob_pcnt"
writeRaster(popindpct, paste0(dir_dat,"/ind/pob-ind.grd"), bylayer=FALSE, format='raster',overwrite=TRUE)
plot(ind)

save.image(paste0(dir_pres,"/raster_processing.RData"))

