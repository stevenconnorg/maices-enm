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
length(todos@data[is.na(todos@data$NomComun),])
length(todos@data[is.na(todos@data$Raza_prima),])

unique(todos@data$Raza_prima)
unique(todos@data$Complejo_r)


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
summary(todos.6)
summary(todos)

# subset races with greater than 20 samples
maices=todos.6



writeOGR(maices,dsn=paste0(dir_maices),layer="todos-maices-cleaned",driver="ESRI Shapefile",overwrite=TRUE)

# maize statistics
library(dplyr)
raza.counts<-maices@data %>% group_by(Raza_prima) %>% summarise(n()) 
raza.counts
unique(maices@data$Raza_prima)
unique(maices@data$Complejo_r)


# make presence/absence matrix
# maices<-readOGR(dsn=paste0(dir_maices,"/todos-maices-cleaned.shp"),layer="todos-maices-cleaned")

# install.packages("letsR")
library(letsR)
xy<-cbind(maices@data$Longitud,maices@data$Latitud)

# extent
# xmin        : -117.625 
# xmax        : -86.20833 
# ymin        : 14.025 
# ymax        : 33.225 
# nrow=2304,ncol=3770

PA<-letsR::lets.presab.points(xy,maices@data$Raza_prima,xmn=-117.625, xmx=-86.20833 , ymn=14.025 , ymx=33.225 , resol = 0.00883333)

plot(PA$Richness_Raster)
pam<-PA$Presence_and_Absence_Matrix
View(pam)
pam[pam == 0 ] <- NA
View(pam)
pa<-data.frame(pam)
pa_df_fn<-paste0(dir_out,"/pa_dataframe.csv")
pam_fn<-paste0(dir_out,"/pam.csv")
con_pa<-file(pa_df_fn,encoding="LATIN1")
con_pam<-file(pam_fn,encoding="LATIN1")

write.csv(pa,file=con_pa)
write.csv(pam,file=con_pam)

save.image(file=paste0(dir_maices,"/clean_maices_obs.RData"))
