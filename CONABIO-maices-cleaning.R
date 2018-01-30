# establish directories
getwd()
root<-"E:/thesis"


dir_dat<-paste0(root,"/01_data")
dir_R<-paste0(root,"/02_R")
dir_out<-paste0(root,"/03_output")
dir_figs<-paste0(root,"/04_figs")
dir_lit<-paste0(root,"/05_lit")
dir_comp<-paste0(root,"/06_comp")
dir_presentations<-paste0(root,"/07_pres")

dir_maices<-paste0(dir_dat,"/maices")
dir_ind<-paste0(dir_dat,"/ind")

# new directories for biomod
dir_bm<-paste0(dir_R,"/01_biomod")
# dir_bmz<-paste0(dir_R,"/02_biomodez")

dir_clim<-paste0(dir_dat,"/clim")
dir_pres<-paste0(dir_clim,"/present")
dir_fut<-paste0(dir_clim,"/future")

dir_p.mosaics<-paste0(dir_pres,"/2.0/")
dir_f.mosaics<-paste0(dir_fut,"/1.4/")

dir_stacks<-paste0(dir_clim,"/stacks/")

folders<-as.list(ls(),pattern="dir")

# tk work on function to create these directories
for (i in folders)  { 
  folder<-get(i)
  dir.create(folder,recursive=TRUE) 
} 

###################################################
# Maize Observation Data Downloading and Cleaning # 
###################################################
dir.create(paste0(dir_maices),recursive=T)
zipurl<-"http://www.conabio.gob.mx/informacion/gis/maps/geo/maizngw.zip"
outdir<-dir_maices
replacement<-"testing"

downloadShapefile <- function(zipurl,outdir,replacement){
  # zipurl as string
  # dir_maices
  # replacement as string
  ext<-file_ext(zipurl)
  filename<-basename(zipurl)
  name<<-gsub(paste0(".",ext),"",filename)
  download.file(zipurl,destfile = paste0(outdir,"/",filename),method = "wget")

  outDir<-normalizePath(outdir,winslash = "\\") # outDir<-paste0(dir_maices) # Define the folder where the zip file should be unzipped to , normalized
  zipF<-paste0(outDir,"/",filename) # lets you choose a file and save its file path in R (at least for windows)
  unzip(zipF,exdir=outDir)  # unzip your file 
  
  files <- list.files(outDir,pattern = name,full.names = F) 
  
  # https://stackoverflow.com/questions/25673643/rename-multiple-files-in-a-folder-using-r

  sapply(files,FUN=function(eachPath){ 
    file.rename(from=eachPath,to= sub(pattern=name, paste0(replacement),eachPath))
  })

}
replacement <- "todos-maices"
downloadShapefile(zipurl,dir_maices,replacement)


library(rgdal)
todos<-readOGR(paste0(dir_maices,"/",replacement,".shp"),layer=replacement,use_iconv=TRUE) 



##############################################
# make sampling bias raster layer 

library(raster) # spatial data manipulation
library(MASS) # for 2D kernel density function
library(magrittr) # for piping functionality, i.e., %>%
library(maptools)
library(rgdal)



setwd(dir_maices)


maiz<-readOGR(dsn="todos-maices.shp",layer = "todos-maices")
todos<-maiz
#for(i in todos@data$i){
#  i<-iconv(todos@data$i, to="LATIN1")
#}
#iconv(todos@data, from="UTF-8", to="LATIN1")

todos@data$Raza_prima<-iconv(todos@data$Raza_prima, from="UTF-8", to="LATIN1")
todos@data$NomComun<-iconv(todos@data$Raza_prima, from="UTF-8", to="LATIN1")
todos@data$Complejo_r<-iconv(todos@data$Complejo_r, from="UTF-8", to="LATIN1")

# inspect data
View(todos@data) 
colnames(todos@data) 
str(todos@data)
head(todos@data$Nom_ent)
head(todos@data$Nom_mun)
head(todos@data$NomComun)
head(todos@data$NomComun)
length(todos@data[is.na(todos@data$NomComun),])
length(todos@data[is.na(todos@data$Raza_prima),])

unique(todos@data$Raza_prima)
unique(todos@data$Complejo_r)
unique(todos@data$AnioColect)


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
todos.6<-todos.5 #[todos.5@data$AnioColect >1969 & todos.5@data$AnioColect < 2001,]# remove period 1 collection

summary(todos.6)
summary(todos)

unique(todos.6@data$Raza_prima)

library(dplyr)

todos.6@data%>%
  group_by(Raza_prima) %>%
  summarise(n = n())%>%
filter(n > 15)

# subset races with greater than 20 samples
maices=todos.6



writeOGR(maices,dsn=paste0(dir_maices),layer="todos-maices-cleaned",driver="ESRI Shapefile",overwrite=TRUE)

# develop KDE sampling bias surface
# following http://harvardforest.fas.harvard.edu/sites/harvardforest.fas.harvard.edu/files/data/p14/hf147/hf147-18-neAnts-4HF.R
bias <- cellFromXY(r, maices)
cells <- unique(sort(bias))
kernelXY <- xyFromCell(r, cells)
samps <- as.numeric(table(bias))

# code to make KDE raster
KDEsur <- sm.density(kernelXY, weights=samps, display="none", ngrid=812, 
                     ylim=c(14.025,33.225), xlim=c(-117.625,-86.20833), nbins=0)
KDErast=SpatialPoints(expand.grid(x=KDEsur$eval.points[,1], y=KDEsur$eval.points[,2]))
KDErast = SpatialPixelsDataFrame(KDErast, data.frame(kde = array(KDEsur$estimate, 
                                                                 length(KDEsur$estimate))))
KDErast <- raster(KDErast)
KDErast <- resample(KDErast, r)
KDErast <- KDErast*r
KDEpts <- rasterToPoints(KDErast)

writeRaster(KDErast,filename=paste0(dir_stacks,"/bias.grd"),overwrite=TRUE)



# maize statistics
library(dplyr)
raza.counts<-maices@data %>% group_by(Raza_prima) %>% summarise(n()) 
raza.counts
unique(maices@data$Raza_prima)
unique(maices@data$Complejo_r)


# make presence/absence matrix
library(rgdal)
maices<-readOGR(dsn=paste0(dir_maices,"/todos-maices-cleaned.shp"),layer="todos-maices-cleaned")
# project maize observations to equal area projection

# proj4string for Mexico Albers Equal Area Conic

# http://spatialreference.org/ref/sr-org/38/



equalarea<-CRS("+proj=aea +lat_1=14.5 +lat_2=32.5 +lat_0=24 +lon_0=-105 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
latlong<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

maices_EA<-spTransform(maices,equalarea)




# add new coordinates to dataframe

maices_EA@data[,c("X_Albers","Y_Albers")]<-coordinates(maices_EA)
maices@data[,c("X_Long","Y_Lat")]<-coordinates(maices)


# extract coordinates as xy object for lets.presab.points

xy_EA<-cbind(maices_EA@data$X_Albers,maices_EA@data$Y_Albers)

xy<-cbind(maices@data$X_Long,maices@data$Y_Lat)

# use template wgs raster for environmental predictors to create presence/absence matrix

r<-raster(nrow=2304,ncol=3770)
extent(r)<-c(-117.625,-86.20833,14.025,33.225)
r
crs(r) <- latlong



# project raster to equal area

rSp<-projectRaster(r,crs=equalarea)

# get extent of object
rSp_extent<-extent(rSp) #(in equalarea projection)
r_extent<-extent(r) #(in equalarea projection)

# xmin        : -1368813
# xmax        : 2031507  
# ymin        : -1117875
# ymax        : 1148009

# use template raster extent as domain for presence/absence matrix
PA_EA<-letsR::lets.presab.points(xy_EA,maices_EA@data$Raza_prima,xmn= rSp_extent@xmin, xmx= rSp_extent@xmax   , ymn=rSp_extent@ymin , ymx= rSp_extent@xmax, resol = res(c(rSp)),crs =equalarea,remove.cells=FALSE,remove.sp=FALSE)
#PA<-letsR::lets.presab.points(xy,maices@data$Raza_prima,xmn= r_extent@xmin, xmx= r_extent@xmax   , ymn=r_extent@ymin , ymx= r_extent@xmax, resol = res(c(r)),crs =latlong,remove.cells=FALSE,remove.sp=FALSE)


# extract matrix from presence/absence object
pam_EA<-PA_EA$Presence_and_Absence_Matrix
#pam<-PA$Presence_and_Absence_Matrix

#View(pam_EA)



# set zero values to NA (not true absences)
pam_EA[pam_EA == 0 ] <- NA
#View(pam_EA)
pa_EA<-data.frame(pam_EA)

#pam[pam == 0 ] <- NA
#View(pam)
#pa<-data.frame(pam)


# # set encoding for file to preserve Spanish encoding
pa_df_fn_EA<-paste0(dir_out,"/pa_dataframe_EA.csv")
pam_fn_EA<-paste0(dir_out,"/pam_EA.csv")

con_pa_EA<-file(pa_df_fn_EA,encoding="LATIN1")
con_pam_EA<-file(pam_fn_EA,encoding="LATIN1")

#pa_df_fn<-paste0(dir_out,"/pa_dataframe.csv")
#pam_fn<-paste0(dir_out,"/pam.csv")

#con_pa<-file(pa_df_fn,encoding="LATIN1")
#con_pam<-file(pam_fn,encoding="LATIN1")


# # write csvs and save images

write.csv(pa_EA,file=con_pa_EA)
write.csv(pam_EA,file=con_pam_EA)

#write.csv(pa,file=con_pa)
#write.csv(pam,file=con_pam)


save.image(file=paste0(dir_maices,"/clean_maices_obs.RData"))

