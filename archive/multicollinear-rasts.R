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
