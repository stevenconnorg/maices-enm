
root<-"/gpfs/home/scg67/thesis"
#root<-"E:/thesis"


setwd(root)

# # #################################################################
# # # NICHE OVERLAP
# # #################################################################
# 
dir_out<-paste0(root,"/03_output")
dir_wmean<-paste0(dir_out,"/wmean-consensus")
dir_binary<-paste0(dir_out,"/binaries")
dir_rangeras<-paste0(dir_out,"/range-size/rasters")
dir_range<-paste0(dir_out,"/range-size")
dir_rangedf<-paste0(dir_out,"/range-size/dfs")
dir_rangeR<-paste0(dir_out,"/range-size/RData")
dir_div<-paste0(dir_out,"/diversity")


library(ENMeval)
m.args<-make.args(make.args(RMvalues = seq(0.5, 5, 0.5), 
                            fc = c("L", "LQ", "H"), 
                            labels = TRUE))
setwd(dir_wmean)
files<-Sys.glob("*current*.grd")
grds<-paste0(getwd(),"/",files)
wmean_current_projs<-stack(grds)

current.overlapI<-calc.niche.overlap(wmean_current_projs, stat = "I", maxent.args)
rownames(current.overlapI)<-sp.n
colnames(current.overlapI)<-sp.n

library(corrplot)
I.mat<-current.overlapI
save(I.mat,file=paste0(dir_div,"/I-stat-current.RData"))
load(paste0(dir_div,"/I-stat-current.RData"))
I.mat[is.na(I.mat)]<-0

png(paste0(dir_div,"/niche_overlap.png"), width =10, height = 13, units = 'in', res = 600)
current.overlapI.corrplot<-corrplot(I.mat, mar = c(0,0,5, 1),hclust.method ="ward.D2", tl.pos = "ld", tl.cex =0.8, tl.col = 'black', method = "square", order="hclust",type = "lower")
dev.off()


# # #################################################################
# # # CALCULATE NICHE SIZE
# # #################################################################
# 


dir_wmean<-paste0(dir_out,"/wmean-consensus")
dir_binary<-paste0(dir_out,"/binaries")
dir_rangeras<-paste0(dir_out,"/range-size/rasters")
dir_range<-paste0(dir_out,"/range-size")dir_rangedf<-paste0(dir_out,"/range-size/dfs")
dir_rangeR<-paste0(dir_out,"/range-size/RData")
dir_div<-paste0(dir_out,"/diversity")


library(ENMeval)

setwd(dir_binary)
projects<-c("proj_current","proj_rcp85_50","proj_rcp85_70")
niche_size<-data.frame(Landrace = as.character(),Niche_Size=as.numeric(),Time=as.character())

for (sp in sp.n){
  for (p in projects[1]){
    bin<-raster(paste0("binary_750_wmean_avg_",sp,"_",p,".grd"))
    p<-sub("proj_","",p)
    bin[is.na(bin)] <- 0
    bin
    dat<-sum(bin[] > 0 ) * res(bin)[1]^2
    dat<-dat/res(bin)[1]^2
    niche_size<-rbind(niche_size,data.frame(sp,dat,p))
  }
}

# convert to square kilometers
niche_size[niche_size=="pro_current"]<-"1970-2000"
niche_size[niche_size=="proj_rcp85_50"]<-"2041-2060"
niche_size[niche_size=="proj_rcp85_70"]<-"2061-280"

write.csv(niche_size,paste0(dir_range,"/niche_size.csv"))

#################################################################
# merge variable importance csvs
#################################################################

setwd(paste0(dir_out,"/var-importance/"))
emvarimp.files <- list.files(pattern="*EM-mean*")
emvarimp.paths<-paste0(getwd(),"/",emvarimp.files)
var.imp.merge<-data.frame()



### merge variable importance (ensemble) csvs into dataframe
merged.em.varimp <-
  do.call(rbind,
          lapply(emvarimp.paths, read.csv))

View(merged.em.varimp)

library(dplyr)
### group variety and variable
dply<-merged.em.varimp %>% group_by_("X","variety")

library(tidyr)

### spread to get columns for each variety
tidy<-dply %>% spread(variety,rowMeans.EMvarimp.)
str(tidy)



### transpose to get variables as columns
tidy.t<-t(tidy)
str(tidy.t)
tidy.df<-as.data.frame(tidy.t)
colnames(tidy.df) <- as.character(unlist(tidy.df[1,]))
tidy.df= tidy.df[-1, ]
View(tidy.df)
str(tidy.df)

tidy.df[] <- lapply(tidy.df, function(x) as.numeric(as.character(x)))
str(tidy.df)

write.csv(tidy.df,paste0(dir_out,"/var-importance/ensemble-merged.csv"))
tidy.df<-read.csv(paste0(dir_out,"/var-importance/ensemble-merged.csv"))
View(tidy.df)

library(corrplot)
var.mat<-as.matrix(tidy.df)

pairs(var.mat)
corr.mat<-cor(t(var.mat))

png(paste0(dir_div,"/var-imp-matrix.png"), width =10, height = 13, units = 'in', res = 600)
corrplot.var<-corrplot(corr.mat,order="hclust",type="upper",diag=FALSE,hclust.method="ward.D2")
dev.off()


library(psych)
cor(t(I.mat), corr.mat)
View(t(I.mat))
View(t(I.mat))

dir_varimp<-paste0(dir_out,"/var-importance")
full.mat<-lowerUpper(corr.mat,upper=I.mat)
ls(pattern="dir")
png(paste0(dir_varimp,"/var-imp-lower-w-niche-overlap-upper.png"), width =10, height = 13, units = 'in', res = 600)
corrplot(full.mat, addrect =5,order="hclust",diag=FALSE,hclust.method="ward.D2",  tl.cex =0.8, tl.col = 'black', method = "square")
dev.off()



# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
head(tidy.df)
rownames(tidy.df)<-tidy.df[,1]
tidy.df<-tidy.df[,-1]
View(tidy.df)
bootvar.df<-t(tidy.df)
full.mat[is.na(full.mat)]<-0

plot(fit2) 

#fit.corr <- pvclust(bootvar.df, method.hclust="ward.D2",method.dist="correlation", nboot=1000,parallel=TRUE)

fit.ward2 <- pvclust(bootvar.df, method.hclust="ward.D2",method.dist="correlation", nboot=1000,parallel=TRUE)
#fit.ward2.full <- pvclust(full.mat, method.hclust="ward.D2",method.dist="correlation", nboot=1000,parallel=TRUE)

?pvrect
png(paste0(dir_div,"/var-imp-cluster-dendrogram.png"), width =10, height = 10, units = 'in', res = 600)

plot(fit.ward2) 

#pvrect(fit.ward2 ,alpha=.985,type="geq")
dev.off()

groups<-cutree(fit.ward2$hclust,3)
sort(groups)
groups<-cutree(fit.ward2$hclust,4)
groups<-cutree(fit.ward2$hclust,5)
#groups<-cutree(fit.ward2$hclust,7)

sort(groups)
class(groups)
save(fit.ward2,file=paste0(dir_div,"/pvclust-var-imp-ward2.RData"))

tidy.df$groups<-groups
colnames(tidy.df)
Anova(tidy.df$ind_pob_pcnt~tidy.df$groups)

tidy.df$groups<-as.factor(tidy.df$groups)
model = lm(tidy.df$BIO15PrecSeas.COV~tidy.df$groups)

aov<-aov(model)
summary(model)
plot(model)


thsd<-TukeyHSD(aov, conf.level = 0.95)
thsd

model2 = lm(tidy.df$EVMPETWarmestQ~tidy.df$groups)

aov2<-aov(model2)
summary(model2)


thsd2<-TukeyHSD(aov2, conf.level = 0.95)
thsd2


varimp.w.group<-tidy.df
write.csv(varimp.w.group,file=paste0(dir_div,"/var-imp-wgrp5.csv"))

means<-varimp.w.group %>% group_by(groups) %>% summarize_all(funs(mean))
sds<-varimp.w.group %>% group_by(groups) %>% summarize_all(funs(sd))

sds
varimpsds<-as.data.frame(sds)
rownames(varimpsds)<-varimpsds[,1]
varimpsds<-varimpsds[,-1]

varimpmeans<-as.data.frame(means)
rownames(varimpmeans)<-varimpmeans[,1]
varimpmeans<-varimpmeans[,-1]
colnames(varimpmeans)
varimpmeans<-varimpmeans[,c(2:12)]
varimpsds<-varimpsds[,c(2:12)]


colnames(varimpmeans)


library(RColorBrewer)
coul = brewer.pal(9, "Pastel2")
colnames(varimpmeans)

png(paste0(dir_out,"/var-importance/env-var-imp-barplot-by-group7.png"), width =10, height = 10, units = 'in', res = 600)
par(mar=c(13, 4.1, 4.1, 2.1))
barCenters<-barplot(as.matrix(varimpmeans, ylim=c(.01,0.4)), main="Variable Importance by Landrace Cluster", ylab="Importance (0-1)",
 	legend = rownames(varimpmeans), las =3, beside=TRUE)
yaxp=c(minY-axis, maxY-axis, Interval)
segments(barCenters, as.matrix(varimpmeans - varimpsds), barCenters,
       as.matrix(varimpmeans  + varimpsds), lwd = 1.5)

arrows(barCenters, as.matrix(varimpmeans - varimpsds[3,12]), barCenters,
       as.matrix(varimpmeans  + varimpsds), lwd = 1.5, angle = 90,
       code = 3, length = 0.05)


dev.off()

ls()


axis(1)
?pvclust
boxplot(tidy.df$ind_pob_pcnt)
tidy.df[tidy.df$ind_pob_pcnt>0.15, c(13,14)]
View(tidy.df)

varimp.means<-tidy.df %>% summarize_all(funs(mean))
t1<-as.data.frame(t(varimp.means))
colnames(t1)[1]<-"Mean"
varimp.sd<-tidy.df %>% summarize_all(funs(sd))
t2<-as.data.frame(t(varimp.sd))
colnames(t2)[1]<-"S.D."
varimp.max<-tidy.df %>% summarize_all(funs(max))
t3<-as.data.frame(t(varimp.max))
colnames(t3)[1]<-"Max"

t1$Variable<-rownames(t1)
t2$Variable<-rownames(t2)
t3$Variable<-rownames(t3)

t5<-merge(t1,t2)
t6<-merge(t5,t3)
View(t6)

write.csv(t3,paste0(dir_out,"/var-importance/mean-sd-max.csv"))

### ECOLOGICAL NICHE OVERLAP STATISTICS
current.overlapI<-calc.niche.overlap(wmean_current_projs, stat = "I", maxent.args)

as.vector(tidy.df[1,])
library(spaa)

var.imp.mat<-tidy.df[,-ncol(tidy.df)]
rm(output)




# c("pianka",  "schoener","petraitis","czech","morisita", "levins")

niche.overlap.mat<-function(input,method) { 
		output<-matrix(, nrow = nrow(input), ncol = nrow(input))
		rownames(output)<-rownames(input)
		colnames(output)<-rownames(input)
		for (i in 1:nrow(input)){
			for (n in 1:nrow(input)){
				stat<-niche.overlap.pair(as.vector(input[i,]),as.vector(input[n,]), method = method)
				output[i,n]<-stat

				}
			}
	return(output)
} 

levins<-niche.overlap.mat(var.imp.mat,method="levins")
schoener<-niche.overlap.mat(input=var.imp.mat,method="schoener")
pianka<-niche.overlap.mat(input=var.imp.mat,method="pianka")
petraitis<-niche.overlap.mat(input=var.imp.mat,method="petraitis")
czech<-niche.overlap.mat(input=var.imp.mat,method="czech")
morisita<-niche.overlap.mat(input=var.imp.mat,method="morisita")

l_dist<-as.matrix(dist(levins))
p_dist<-as.matrix(dist(pianka))

corrplot(cor(p_dist),type="lower", method = "square")

econiche.mat<-lowerUpper(p_dist,upper=t(l_dist))

corrplot(econiche.mat, corr=,method = "square")

save(list=c("pianka",  "schoener","petraitis","czech","morisita", "levins"),file=paste0(dir_div,"/niche_overlap_statistics.RData"))
#########################
setwd(paste0(dir_out,"/var-importance/"))
varimp.files <- list.files(pattern="*model-var-imp*")
varimp.paths<-paste0(getwd(),"/",varimp.files)
var.imp.merge<-data.frame()


### merge individual variable importance csvs into dataframe
       merged.varimp <-
         do.call(rbind,
                 lapply(varimp.paths, read.csv))
       library(dplyr)
library(reshape)
var.melt<-merged.varimp %>% melt(id=c("X","variety"))
View(var.melt)

var.melt.s<-var.melt %>% spread(key = X,value=value)

View(var.melt.s)
avg.imp<-var.melt.s %>% group_by(variety) %>% summarize_all(funs(mean))
avg.imp<-as.data.frame(avg.imp)
rownames(avg.imp)<-avg.imp[,1]
avg.imp<-avg.imp[,c(-1,-2)]
str(avg.imp)

      # 
      # 
      # ### transpose to get variables as columns
      # View(dply.2)
      # colnames(dply.2)[1]<-"Variable"
      # rownames(dply.2)<-dply.2$X
      # str(dply.2)
      # dply.2[,2]<-as.numeric(as.character(dply.2[,2]))
      # dply.2[,3]<-as.numeric(as.character(dply.2[,3]))
      # dply.2[,4]<-as.numeric(as.character(dply.2[,4]))
      # dply.2[,5]<-as.numeric(as.character(dply.2[,5]))
      # dply.2[,6]<-as.numeric(as.character(dply.2[,6]))
      # dply.2[,7]<-as.numeric(as.character(dply.2[,7]))
      # dply.2[,8]<-as.numeric(as.character(dply.2[,8]))
      # 
      # View(dply.2)
      # class(dply.2)
      # colnames(dply.2)
      # 
      # dply.3<-t(dply.2)
      # View(dply.3)
      # model = lm(Mean ~ Variable * variety,
      #            data=dply.2)
      # 
      # model 
      # # 
      # library(car)
      # 
      # aov<-aov(model, Type="II")
      # aov
      # Anova<-Anova(model, Type="II",
      #              white.adjust=TRUE)
      # Anova
      # 
      # 
      # thsd<-TukeyHSD(aov,ordered = TRUE, conf.level = 0.95)
      # thsd
      # aov()
      # colnames(dply.2)
      # summary(tidy2.df)
      # plot(tidy2.df,type="l")
      # 
      # var.pca.d


# 
# # #################################################################
# # # merge model evaluation csvs
# # #################################################################
# 
setwd(paste0(dir_out,"/ensemble-evals/"))
emval.files <- list.files(pattern="*EM_eval*")
emval.paths<-paste0(getwd(),"/",emval.files)

# ### merge variable importance (ensemble) csvs into dataframe
merged.em.evals <-
  do.call(rbind,
          lapply(emval.paths, read.csv))

View(merged.em.evals)
colnames(merged.em.evals)


### group variety and variable
dply<-merged.em.evals %>% group_by(variety,Model.name,Eval.metric) %>% summarize(Testing.data=mean(Testing.data),Sensitivity=mean(Sensitivity),Specificity=mean(Specificity))
View(dply)

dply<-as.data.frame(dply)
dply$Model.name<-sub("_MAXENT.Phillips_mergedRun_mergedData","",dply$Model.name)

var.avg<-dply %>% group_by(variety,Eval.metric,Model.name) %>% summarize(Testing.data=mean(Testing.data),Sensitivity=mean(Sensitivity),Specificity=mean(Specificity))
View(var.avg)

# model scores tend to be similar for each wmean ensemble
# average by evaluation metric
var.avg2<-dply %>% group_by(variety,Eval.metric) %>% summarize(Testing.data=mean(Testing.data),Sensitivity=mean(Sensitivity),Specificity=mean(Specificity))
View(var.avg2)

eval.metrics.by.var<-as.data.frame(var.avg2)
colnames(eval.metrics.by.var)
head(eval.metrics.by.var)
Kappa_by_sp<-subset(eval.metrics.by.var,Eval.metric %in% c("KAPPA"))
sp.n<-unique(Kappa_by_sp$variety)
tau<-read.csv(paste0(dir_maices,"/tau.csv"))
head(tau)
tau$Raza_prima<-sub(" ",".",tau$Raza_prima)
tau$Raza_prima<-sub(" ",".",tau$Raza_prima)
tau$Raza_prima<-sub(" ",".",tau$Raza_prima)
tau$Raza_prima<-sub(" ",".",tau$Raza_prima)
tau$Raza_prima

sp_tau<-tau[tau$Raza_prima %in% sp.n,]

tau_kappa<-merge(sp_tau,Kappa_by_sp,by.x="Raza_prima",by.y="variety")
nichesize_kappa<-merge(Kappa_by_sp,niche_size,by.x="variety",by.y="sp")
head(nichesize_kappa)
str(nichesize_kappa)
View(tau_kappa)
colnames(tau_kappa)
cor(nichesize_kappa$Testing.data,nichesize_kappa$dat)

niche_size
write.csv(var.avg2,file=paste0(dir_out,"/ensemble-evals/avg_score_by_metric_by_variety_with_ss.csv"))

# get mean sensitivity and specificity by variety
meanss<-var.avg2 %>% group_by(variety) %>% summarize(mean(Sensitivity),mean(Specificity))
meanss<-as.data.frame(meanss)
mean(meanss$`mean(Sensitivity)`)
#  78.5941
mean(meanss$`mean(Specificity)`)
#  91.15314
write.csv(meanss,file=paste0(dir_out,"/ensemble-evals/mean-spec-sens-by-variety.csv"))

spread(var.avg,Eval.metric,Testing.data)

### spread to get columns for each variety
colnames(var.avg)
var.avg2<-var.avg[,1:3]
tidy<-var.avg %>% group_by(variety,Eval.metric) %>% summarize(mean(Sensitivity),mean(Specificity)) %>% spread(key=c(Eval.metric),value=Testing.data)
View(tidy)
auc<-read.csv("E:\\thesis\\03_output\\ensemble-evals\\mean_auc_by_variety.csv")
tidy<-as.data.frame(tidy)
tidyauc<-merge(tidy,auc,by="variety")
mean(tidyauc$Mean_AUC)
# mean AUC 0.9392401
write.csv(tidyauc,file=paste0(dir_out,"/ensemble-evals/avg_score_by_variety_spread_by_metric-w-auc.csv"))
# 
# # #################################################################
# # # merge range change dataframe
# # #################################################################
# # 

setwd(paste0(dir_out,"/range-size/dfs/binary_750_wmean_avg_"))
chg.files <- list.files(pattern="*")
chg.paths<-paste0(getwd(),"/",chg.files)

# ### merge variable importance (ensemble) csvs into dataframe
merged.rangesize <-
  do.call(rbind,
          lapply(chg.paths, read.csv))
library(tidyr)
library(dplyr)
colnames(merged.rangesize)
View(merged.rangesize)
colnames(merged.rangesize)

rangesorted <- merged.rangesize[order(merged.rangesize$time, merged.rangesize$species),]### group variety and variable


View(rangesorted)
colnames(rangesorted)


percchang<-rangesorted[,c("species","time","SpeciesRangeChange")]
percchang %>% select(SpeciesRangeChange,species,time)%>% spread(key=c(time),value=SpeciesRangeChange, convert=TRUE)
View(percchang)

colnames(percchang)[1]<-"variety"
colnames(percchang)[2]<-"Range_Change_%_(2041-2060)"
colnames(percchang)[3]<-"Range_Change_%_(2061-2070)"

write.csv(percchang,file=paste0(dir_range,"/range-change-by-period.csv"))
View(percchang2050)

# 
# 
# #################################################################
# # diversity zonal stats
# #################################################################
# #import required libraries
projects<-c("proj_current","proj_rcp85_50","proj_rcp85_70")
var.impgrp<-read.csv(paste0(dir_div,"/var-imp-wgrp5.csv"))
colnames(var.impgrp)
groups<-var.impgrp$groups
binaries<-Sys.glob(paste0(dir_binary,"/binary_750_wmean*proj_current.grd"))
temp<-raster(binaries[1])

maskmex<-readOGR(paste0(dir_dat,"/mexico-outline.shp"))
 maskmex<-spTransform(maskmex, crs(temp))

bin_stack<-stack()

		for (p in projects){
			binaries<-Sys.glob(paste0(dir_binary,"/binary_750_wmean*",p,".grd"))
				for (bin in binaries){
					bin<-raster(bin)
					bin<-resample(bin,temp)
					bin_stack<-stack(bin_stack,bin)}
					#writeRaster(bin_stack,file=paste0(dir_div,"/binary_stack_",p,".tif"),format="GTiff",overwrite=TRUE)
				      div_by_group<-stackApply(bin_stack,indices=groups,fun="sum")
					writeRaster(div_by_group,file=paste0(dir_div,"/div_by_group_",p,".grd"))

bin_stack<-stack()

}


		for (p in projects){
				bin_stack<-stack(paste0(dir_div,"/binary_stack_",p,".grd"))
				div_by_group<-stackApply(bin_stack,indices=groups,fun="sum")

				names(div_by_group)<-paste0("Group_", seq(1,nlayers(div_by_group),1))
				div_by_group@title<-p

				div_by_group<-mask(div_by_group,maskmex,inverse=FALSE)


				writeRaster(div_by_group,file=paste0(dir_div,"/div_by_group_",p,".grd"),
				overwrite=TRUE)
}

}
cdiv_stack<-stack(paste0(dir_div,"/div_by_group_",projects[1],".grd"))


## get diversity change by group
for (p in projects[2:3]){
  div_stack<-stack(paste0(dir_div,"/div_by_group_",p,".grd"))
  for (i in 1:nlayers(div_stack)){
    div_stack[[i]] = div_stack[[i]] -  cdiv_stack[[i]] 
  }
  writeRaster(div_stack,file=paste0(dir_div,"/div_by_group_chg_",p,".grd"),overwrite=TRUE)
}

##### EVALUATE DIVERSITY BY INDIGENOUS PRESENCE #####

# dir_range<-paste0(dir_out,"/range-size")
# #list files (in this case raster TIFFs)

alpha.div <- raster(paste0(dir_div,"/maize_alpha_diversity_current.grd"))
# alpha.div50 <- raster(paste0(dir_div,"/maize_alpha_diversity_rcp85_50.grd"))
# alpha.div70 <- raster(paste0(dir_div,"/maize_alpha_diversity_rcp85_70.grd"))

alpha.divgrp <- stack(paste0(dir_div,"/div_by_group_proj_current.grd"))

alpha.div.grp.masked<-mask(alpha.divgrp,mask=maskmex,inverse=F)
alpha.div.masked<-mask(alpha.div,mask=maskmex,inverse=F)

writeRaster(alpha.div.grp.masked,paste0(dir_div,"/div_by_group_proj_current-masked.grd"))
writeRaster(alpha.div.masked,paste0(dir_div,"/maize_alpha_diversity_current-masked.grd"))



col.l <- colorRampPalette(c('blue', 'green',  'yellow', 'red'))(1000)


names(alpha.div.grp.masked)<-paste0("Group_",seq(1,5,1))
calphadivlevel<-levelplot(alpha.div.grp.masked,col.regions=col.l, names.attr=names(alpha.div.grp.masked),maxpixels=2e5, margin = list(FUN = 'mean'),
                         xlab="X", ylab="Y",main=paste0("Alpha Diversity "))
calphadivlevel

      png(paste0(dir_div,"/alpha-div-grp.png"), width = 10, height = 10, units = 'in', res = 600)
      print(calphadivlevel)
      dev.off()


calphadiv<-levelplot(alpha.div.masked,col.regions=col.l,maxpixels=2e5, margin = list(FUN = 'mean'),
                         xlab="X", ylab="Y",main=paste0("Alpha Diversity 1970-2000"))
calphadiv

      png(paste0(dir_div,"/current-alpha-div.png"), width = 10, height = 10, units = 'in', res = 600)
      print(calphadiv)
      dev.off()


# 
dir_wmean<-paste0(dir_out,"/wmean-consensus")
dir_binary<-paste0(dir_out,"/binaries")
dir_rangeras<-paste0(dir_out,"/range-size/rasters")
dir_range<-paste0(dir_out,"/range-size")
dir_rangedf<-paste0(dir_out,"/range-size/dfs")
dir_rangeR<-paste0(dir_out,"/range-size/RData")
dir_div<-paste0(dir_out,"/diversity")


library(maptools)
library(raster)
library(rgdal)
library(rasterVis)

# read in alpha diversity by group stack

for (p in projects){
divs<-stack(paste0(dir_div,"/div_by_group_",p,".grd"))
cdiv_bg<-mask(stack1 ,maskmex,inverse=F)
}

plot(maskmex)

cdiv_bg.m<-mask(cdiv_bg ,maskmex,inverse=T)
levelplot(cdiv_bg)

# #read-in the municipios polygon shapefile
poly <- readOGR("E:\\thesis\\01_data\\ind\\presindigw\\presindigw.shp",layer="presindigw")
ob <- SpatialPolygons(poly@polygons,proj4string=poly@proj4string)
ob<-spTransform(ob,proj4string(cdiv_bg))
poly<-spTransform(poly,proj4string(cdiv_bg))

# get dataframe
df<-poly@data
df<-as.data.frame(df)
df.by.grp<-df


for (p in projects){
cdiv_bg<-stack(paste0(dir_div,"/div_by_group_",p,".grd"))
names(cdiv_bg)<-paste0("Group_",seq(1,nlayers(cdiv_bg),1))



for (lay in 1:nlayers(cdiv_bg)){

	# #extract raster cell count (sum) within each polygon area (poly)
	library(spatialEco)

	mean.stat <- function(x, na.rm=TRUE) {
	   if (na.rm)
	       x <- x[!is.na(x)]
	  sum(x)/length(x)
	  }

	med.stat <- function(x, na.rm=TRUE) {
	   if (na.rm)
	       x <- x[!is.na(x)]
		median(x)
	  }

	min.stat <- function(x, na.rm=TRUE) {
	   if (na.rm)
	       x <- x[!is.na(x)]
		min(x)
	  }
	# 

	mean.div<-zonal.stats(poly, cdiv_bg[[lay]], stat=mean.stat, plot = FALSE)
	df.by.grp$mean.div<-mean.div
	colnames(df.by.grp)[ncol(df.by.grp)] <- paste0(names(cdiv_bg[[lay]]),"_",p,"_mean_div")

}
}

write.csv(df.by.grp, file = paste0(dir_out,"/alpha-div-mean-by-grp.csv"))


colnames(df.by.grp)
df.by.grp<-df.by.grp[,-c(24:28)]
for (i in c(19:23)){
chg<-((df.by.grp[,i]-df[,i-5]))
df.by.grp$chg<-chg
colname<-colnames(df[i])
newname<-sub("_mean_div","_div_chg",colname)
colnames(df.by.grp)[ncol(df)] <- newname
}

colnames(df.by.grp)[18:22]
str(df)

write.csv(df, file = paste0(dir_out,"/alpha-div-mean-by-grp-2050chg.csv"))
df<-read.csv(file = paste0(dir_out,"/alpha-div-mean-by-grp-2050.csv"))
colnames(df)

#df.50chg<-read.csv(file = paste0(dir_div,"/alpha-div-mean-by-grp-df-2050chg.csv"))
df<-read.csv(file = paste0(dir_div,"/alpha-div-by-muni-df.csv"))



levels(df$Presencia)

View(df)
# 
str(df)
df$Presencia<-factor(df$Presencia)
levels(df$Presencia)
levels(df$marmun)


df$Presencia<-factor(df$Presencia,levels=c("Sin población indí­gena","Población indí­gena dispersa","Población con presencia indí­gena","Población indí­gena"),order=TRUE)
df$Presencia<-factor(df$Presencia,levels=c(levels(df$Presencia)[4],levels(df$Presencia)[3],levels(df$Presencia)[1],levels(df$Presencia)[2]),order=TRUE)
df$marmun<-factor(df$marmun,levels=c("Muy bajo","Bajo","Medio","Alto","Muy alto"),order=TRUE)



colnames(df)

### TEST DIVERSITY STATS BY INDIGENOUS POP PRESENCE/MARGINALITY
model = lm(mean.div~ df$Presencia,
           data=df)
model2 = lm(mean.div~ df$marmun,
           data=df)

library(car)

## try by indigenous presence
bartlett<-bartlett.test(mean.div~ df$Presencia,data=df)
bartlett
fligner<-fligner.test(mean.div~ df$Presencia,data=df)
fligner
leveneTest<-leveneTest(mean.div~ df$Presencia,data=df)
leveneTest


# heteroskedasticity present

df<-na.omit(df)
colnames(df)
oneway<-oneway.test(df$med.div~Presencia, data=df, na.action=na.omit, var.equal=FALSE)


TukeyHSD<-TukeyHSD(aov(med.div~ df$Presencia,
           data=df))
TukeyHSD

## try by marginality level
bartlett<-bartlett.test(med.div~ df$marmun,data=df)
bartlett
fligner<-fligner.test(med.div~ df$marmun,data=df)
fligner
leveneTest<-leveneTest(med.div~ df$marmun,data=df)
leveneTest
# heteroskedasticity also present

capture.output(bartlett,file=paste0(dir_out,"/ind/bartlett_med-div_by_marginality.txt"))
capture.output(fligner,file=paste0(dir_out,"/ind/fligner_med-div_by_marginality.txt"))
capture.output(leveneTest,file=paste0(dir_out,"/ind/leveneTest_med-div_by_marginality.txt"))

colnames(df)
oneway<-oneway.test(med.div~marmun, data=df, na.action=na.omit, var.equal=FALSE)
capture.output(oneway,file=paste0(dir_out,"/ind/onewayunequal_med-div_by_marginality.txt"))

str(oneway)
TukeyHSD<-TukeyHSD(aov(med.div~ df$marmun,
           data=df))


str(TukeyHSD[df$marmun])

write.csv(TukeyHSD$`df$marmun`,paste0(dir_out,"/ind/tukeyhsd_med-div_by_marginality.csv"))

kruskal<-kruskal.test(med.div~marmun, data=df, na.action=na.omit)
kruskalp<-kruskal.test(med.div~Presencia, data=df, na.action=na.omit)

library(pastecs)  # for descriptive statistics
library(pgirmess)
install.packages("pgirmess")
# creating a rank variable
df$med.divRank <- rank(df$med.div)

# getting the descriptive statistics for the groups
by(df$med.divRank, df$Presencia, stat.desc, basic = FALSE)
by(df$med.divRank, df$marmun, stat.desc, basic = FALSE)


# post-hoc test for identifying statistical significant differences between the groups
kruskalmc(med.div~ Presencia, data = df)
kruskalmc(med.div~ marmun, data = df)



capture.output(kruskal,file=paste0(dir_out,"/ind/kruskal_med-div_by_marginality.txt"))

df %>% group_by(marmun) %>% summarize(mean(med.div))


# https://rpubs.com/aaronsc32/games-howell-test
games.howell <- function(grp, obs) {
  
  #Create combinations
  combs <- combn(unique(grp), 2)
  
  # Statistics that will be used throughout the calculations:
  # n = sample size of each group
  # groups = number of groups in data
  # Mean = means of each group sample
  # std = variance of each group sample
  n <- tapply(obs, grp, length)
  groups <- length(tapply(obs, grp, length))
  Mean <- tapply(obs, grp, mean)
  std <- tapply(obs, grp, var)
  
  statistics <- lapply(1:ncol(combs), function(x) {
    
    mean.diff <- Mean[combs[2,x]] - Mean[combs[1,x]]
    
    #t-values
    t <- abs(Mean[combs[1,x]] - Mean[combs[2,x]]) / sqrt((std[combs[1,x]] / n[combs[1,x]]) + (std[combs[2,x]] / n[combs[2,x]]))
    
    # Degrees of Freedom
    df <- (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]])^2 / # Numerator Degrees of Freedom
            ((std[combs[1,x]] / n[combs[1,x]])^2 / (n[combs[1,x]] - 1) + # Part 1 of Denominator Degrees of Freedom 
              (std[combs[2,x]] / n[combs[2,x]])^2 / (n[combs[2,x]] - 1)) # Part 2 of Denominator Degrees of Freedom
    
    #p-values
    p <- ptukey(t * sqrt(2), groups, df, lower.tail = FALSE)
    
    # Sigma standard error
    se <- sqrt(0.5 * (std[combs[1,x]] / n[combs[1,x]] + std[combs[2,x]] / n[combs[2,x]]))
          
    # Upper Confidence Limit
    upper.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff + qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Lower Confidence Limit
    lower.conf <- lapply(1:ncol(combs), function(x) {
      mean.diff - qtukey(p = 0.95, nmeans = groups, df = df) * se
    })[[1]]
    
    # Group Combinations
    grp.comb <- paste(combs[1,x], ':', combs[2,x])
    
    # Collect all statistics into list
    stats <- list(grp.comb, mean.diff, se, t, df, p, upper.conf, lower.conf)
  })
  
  # Unlist statistics collected earlier
  stats.unlisted <- lapply(statistics, function(x) {
    unlist(x)
  })
  
  # Create dataframe from flattened list
  results <- data.frame(matrix(unlist(stats.unlisted), nrow = length(stats.unlisted), byrow=TRUE))
  
  # Select columns set as factors that should be numeric and change with as.numeric
  results[c(2, 3:ncol(results))] <- round(as.numeric(as.matrix(results[c(2, 3:ncol(results))])), digits = 3)
  
  # Rename data frame columns
  colnames(results) <- c('groups', 'Mean Difference', 'Standard Error', 't', 'df', 'p', 'upper limit', 'lower limit')

  return(results)
}

df<-na.omit(df)

gh<-games.howell(df$marmun, df$med.div)

write.csv(gh,paste0(dir_out,"/ind/games_howell_med-div_by_marginality.csv"))


colnames(df)

gh<-games.howell(df$Presencia, df$med.div)
write.csv(gh,paste0(dir_out,"/ind/games_howell_med-div_by_indpres.csv"))



save(list=c("gh","bartlett","fligner","leveneTest","oneway","TukeyHSD","df"),file=paste0(dir_div,"/one-way-by-ind-presence.RData"))


bartlett<-bartlett.test(df$med.div, df$Presencia)
fligner<-fligner.test(df$med.div, df$Presencia)
leveneTest<-leveneTest(df$med.div, df$Presencia)
# heteroskedasticity present



capture.output(bartlett,file=paste0(dir_out,"/ind/bartlett_med-div_by_indpres.txt"))
capture.output(fligner,file=paste0(dir_out,"/ind/fligner_med-div_by_indpres.txt"))
capture.output(leveneTest,file=paste0(dir_out,"/ind/leveneTest_med-div_by_indpres.txt"))

colnames(df)
oneway<-oneway.test(med.div~Presencia, data=df, na.action=na.omit, var.equal=FALSE)
capture.output(oneway,file=paste0(dir_out,"/ind/onewayunequal_med-div_by_indpres.txt"))

str(oneway)
TukeyHSD<-TukeyHSD(aov(med.div~ df$Presencia,
           data=df))


write.csv(TukeyHSD$`df$Presencia`,paste0(dir_out,"/ind/tukeyhsd_med-div_by_indpres.csv"))

oneway<-oneway.test(med.div~Presencia, data=df, na.action=na.omit, var.equal=FALSE)

TukeyHSD<-TukeyHSD(aov(med.div~ Presencia, data = df))
gh<-games.howell(df$Presencia, df$med.div)


kruskal<-kruskal.test(med.div~Presencia, data=df, na.action=na.omit)
capture.output(kruskal,file=paste0(dir_out,"/ind/kruskal_med-div_by_indpres.txt"))
bf.test(med.div~Presencia, df, alpha = 0.05, na.rm = TRUE, verbose = TRUE)
bf.test(med.div~marmun, df, alpha = 0.05, na.rm = TRUE, verbose = TRUE)

colnames(df)
df$lg1<-log(df$Group_1_mean_div+1)
df$lg2<-log(df$Group_2_mean_div+1)
df$lg3<-log(df$Group_3_mean_div+1)
df$lg4<-log(df$Group_4_mean_div+1)
df$lg5<-log(df$Group_5_mean_div+1)

qqnorm(df$lg1)
qqnorm(df$lg2)
qqnorm(df$lg3)
qqnorm(df$lg4)
qqnorm(df$lg5)


# then build models for each group
shapiro.test(df$lg1)
shapiro.test(df$lg2)
shapiro.test(df$lg3)
shapiro.test(df$lg4)
shapiro.test(df$lg5)
hist(df$lg1)

# then build models for each group
bartlett<-bartlett.test(df$lg1, df$Presencia)
bartlett<-bartlett.test(df$lg2, df$Presencia)
bartlett<-bartlett.test(df$lg3, df$Presencia)
bartlett<-bartlett.test(df$lg4, df$Presencia)
bartlett<-bartlett.test(df$lg5, df$Presencia)

# heteroskedasticity present
hist(df$Group_1_mean_div)
hist(df$Group_2_mean_div)
hist(df$Group_3_mean_div)
hist(df$Group_4_mean_div)
hist(df$Group_5_mean_div)
View(df
df$indpresmarg <- with(df, interaction(marmun,  Presencia), drop = TRUE )
head(df$indpresmarg)
kruskal1<-kruskal.test(lg1~marmun, data=df, na.action=na.omit)
kruskal2<-kruskal.test(lg2~marmun, data=df, na.action=na.omit)
kruskal3<-kruskal.test(lg3~marmun, data=df, na.action=na.omit)
kruskal4<-kruskal.test(lg4~marmun, data=df, na.action=na.omit)
kruskal5<-kruskal.test(lg5~marmun, data=df, na.action=na.omit)

kruskal1p<-kruskal.test(Group_1_mean_div~Presencia, data=df, na.action=na.omit)
kruskal2p<-kruskal.test(Group_2_mean_div~Presencia, data=df, na.action=na.omit)
kruskal3p<-kruskal.test(Group_3_mean_div~Presencia, data=df, na.action=na.omit)
kruskal4p<-kruskal.test(Group_4_mean_div~Presencia, data=df, na.action=na.omit)
kruskal5p<-kruskal.test(Group_5_mean_div~Presencia, data=df, na.action=na.omit)


kruskal1
kruskal2
kruskal3
kruskal4
kruskal5

install.packages("FSA")
library(FSA)

PT1 = dunnTest(lg1~ marmun,
              data=df,
              method="bh")
PT2 = dunnTest(lg2~ marmun,
              data=df,
              method="bh")
PT3 = dunnTest(lg3~ marmun,
              data=df,
              method="bh")
PT4 = dunnTest(lg4~ marmun,
              data=df,
              method="bh")
PT5 = dunnTest(lg5~ marmun,
              data=df,
              method="bh")
colnames(df)
gh1<-games.howell(df$Presencia, df$Group_1_mean_div)
gh2<-games.howell(df$Presencia, df$Group_2_mean_div)
gh3<-games.howell(df$Presencia, df$Group_3_mean_div)
gh4<-games.howell(df$Presencia, df$Group_4_mean_div)
gh5<-games.howell(df$Presencia, df$Group_5_mean_div)

install.packages("onewaytests")
library(onewaytests)
bf.test(Group_1_mean_div ~ Presencia, df, alpha = 0.05, na.rm = TRUE, verbose = TRUE)
bf.test(Group_2_mean_div ~ Presencia, df, alpha = 0.05, na.rm = TRUE, verbose = TRUE)
bf.test(Group_3_mean_div ~ Presencia, df, alpha = 0.05, na.rm = TRUE, verbose = TRUE)
bf.test(Group_4_mean_div ~ Presencia, df, alpha = 0.05, na.rm = TRUE, verbose = TRUE)
bf.test(Group_5_mean_div ~ Presencia, df, alpha = 0.05, na.rm = TRUE, verbose = TRUE)

gh1m<-games.howell(df$marmun, df$Group_1_mean_div)
gh2m<-games.howell(df$marmun, df$Group_2_mean_div)
gh3m<-games.howell(df$marmun, df$Group_3_mean_div)
gh4m<-games.howell(df$marmun, df$Group_4_mean_div)
gh5m<-games.howell(df$marmun, df$Group_5_mean_div)

gh1m
gh2m
gh3m
gh4m
gh5m

# then build models for each group
df$marmuno<-factor(df$marmun, ordered = FALSE )
df$Presenciao<-factor(df$Presencia, ordered = FALSE )
colnames(df)

df$indpcnt<-df$pobindi/df$pobtot


lm1<-lm(Group_1_mean_div~ marmun+Presencia+indpcnt, data = df)
lm2<-lm(Group_2_mean_div~ marmun+Presencia+indpcnt, data = df)
lm3<-lm(Group_3_mean_div~ marmun+Presencia+indpcnt, data = df)
lm4<-lm(Group_4_mean_div~ marmun+Presencia+indpcnt, data = df)
lm5<-lm(Group_5_mean_div~ marmun+Presencia+indpcnt, data = df)

summary(lm1)
summary(lm2)
summary(lm3)
summary(lm4)
summary(lm5)


colnames(df)
lm1<-lm(med.div ~ marmuno+Presenciao+pobindi, data = df)


lm1
plot(lm1)
summary(lm1)
qqnorm(lm1$residual)
plot(lm1$residual)

summary(lm2)
summary(lm3)
summary(lm4)
summary(lm5)


as.matrix(man.dat[,1:4])
RVAideMemoire::mshapiro.test(as.matrix(man.dat[,1:4]))
y2,y3)
fit <- manova(Y ~ A*B)

res.man <- manova(mean.divs ~ df$Presencia, data = df)
summary(res.man, test="Pillai")


View(gh)

# get categories of indigenous population levels
df$indpobpct<- df$pobindi/df$pobtot
df$indpobpct<-as.numeric(df$indpobpct)
df$ind_pct_group <- cut(df$indpobpct, 5)

# str(df)
# nrow(poly@data)
# nrow(df)
# 
# poly2<-poly
# poly2@data<-merge(poly@data,df)
# nrow(poly2@data)
# length(unique(poly2@data$MUN_OFICIA))
# poly2@data$ID<-c(1:nrow(poly2@data))
# colnames(poly2@data)
# 
# 
# colnames(df)
# colnames(df)
# as.numeric(df$mean.div)
# df<-as.data.frame(df)
# View(df)
# nrow(df)
# 
# save.image(file="E:\\thesis\\01_data\\div_aov.RData")
# 
# 
# 
# 
# 
df
library(rgdal)
presind <- readOGR("E:\\thesis\\01_data\\ind\\poinmun10gw\\poinmun10gw.shp",layer="poinmun10gw")
presind <- readOGR("E:\\thesis\\01_data\\ind\\presindigw\\presindigw.shp",layer="presindigw")
presind@data<-merge(presind@data,df,by="MUN_OFICIA")
writeOGR(obj=presind , dsn="E:\\thesis\\01_data\\ind\\diveristy_by_municipio.shp", layer="diveristy_by_municipio", driver="ESRI Shapefile") 

# 
# 
# View(df.2)
# 
# #write to a CSV file
# names(df.2) <- tolower(names(df.2))
# 
# str(df.2)
# names(df.2) <- tolower(names(df.2))
# library(leaps)
# 
# 
# regsubsets.out <-
#     regsubsets(min.div ~ p3t10+ p3i10+ mon10+ bil10+ nei10+  nli10+ neli10+area.x+perimeter.x+presencia + pobindi+marmun ,
#                data = df.2,
#                nbest = 1,       # 1 best model for each number of predictors
#                nvmax = NULL,    # NULL for no limit on number of variables
#                force.in = NULL, force.out = NULL,
#                method = "exhaustive")
# regsubsets.out
# 
# summary.out <- summary(regsubsets.out)
# as.data.frame(summary.out$outmat)
# 
# plot(regsubsets.out, scale = "adjr2", main = "Adjusted R^2")
# library(car)
# layout(matrix(1:2, ncol = 2))
# ## Adjusted R2
# res.legend <-
#     subsets(regsubsets.out, statistic="adjr2", legend = FALSE, min.size = 5, main = "Adjusted R^2")
# ## Mallow Cp
# res.legend <-
#     subsets(regsubsets.out, statistic="cp", legend = FALSE, min.size = 5, main = "Mallow Cp")
# 
# which.max(summary.out$adjr2)
# 
# best.model <- lm(mea.div ~ p3t10+ p3i10+  nei10+  neli10+area.x+perimeter.x+presencia*marmun+ pobindi , data = df.2)
# summary(best.model)
# 
# 
# summary.out$which[11,]
# abline(a = 1, b = 1, lty = 2)
# library(corrgram)
# library("ggpubr")
# ggscatter(df, x = "index_1", y = "pobpct", 
#           add = "reg.line", conf.int = TRUE, 
#           cor.coef = TRUE, cor.method = "pearson",
#           xlab = "Maize Alpha-diversity", ylab = "Indigenous Population (%)")
# 
# 
# # Shapiro-Wilk normality tests 
# shapiro.test(df$index_1) # => p = 0.1229
# 
# shapiro.test(df$pobpct) # => p = 0.09
# 
# library("ggpubr")
# # mpg
# ggqqplot(df$index_1, ylab = "Alpha Diversity")
# # wt
# ggqqplot(df$pobpct, ylab = "Ind. Pob. (%)")
# 
# cor.pearson <-cor.test(df$index_1,df$pobpct,  method = c("pearson"))
# cor.pearson
# cor.kendall <-cor.test(df$index_1,df$pobpct,  method = "kendall")
# cor.kendall
# cor.spearman <-cor.test(df$index_1,df$pobpct,  method = "spearman")
# cor.spearman
# 
# 
# 
# # Figner-Killeen Test of Homogeneity of Variances
# colnames(df)
# str(df)
# 
# head(df$group)
# fligner.test(index_1~group, data=df)
# bartlett.test(index_1~group, data = df)
# kruskal.test(index_1~group, data = df)
# View(df$index_1)
# 
# 
# 
# plot(alpha.div)
# pobindras<-raster("E:\\thesis\\01_data\\ind\\pob-ind.grd")
# library(spatialEco)
# pobindras<-projectRaster(pobindras,alpha.div)
# pobindras<-resample(pobindras,alpha.div,method="bilinear")
# 
# rascorr3<-rasterCorrelation(pobindras, alpha.div, s = 3, type = "pearson", file.name = NULL)
# rascorr3<-rascorr
# plot(rascorr3)
# rascorr5<-rasterCorrelation(pobindras, alpha.div, s = 5, type = "pearson", file.name = NULL)
# rascorr7<-rasterCorrelation(pobindras, alpha.div, s = 7, type = "pearson", file.name = NULL)
# rascorr25<-rasterCorrelation(pobindras, alpha.div, s = 25, type = "pearson", file.name = NULL)
# 
# 
# library(raster)
# library(SDMTools)
# library(adehabitat)
# 
# rAsc <- asc.from.raster(pobindras) # Function from SDMTools to convert to asc format
# pobindrasspdf <- asc2spixdf(rAsc)
# 
# rAsc <- asc.from.raster(alpha.div) # Function from SDMTools to convert to asc format
# alpha.divspdf<- asc2spixdf(rAsc)
# 
# class(alpha.div)
# modttest<-spatialEco::raster.modifed.ttest(pobindrasspdf , alpha.divspdf, d = "AUTO",
#   sub.sample = TRUE, type = "hexagon", p = 0.1)
# 
# modttest






#mean.div<-zonal.stats(poly, alpha.div[[lay]], stat=mean.stat, trace = TRUE, plot = FALSE)
#med.div<-zonal.stats(poly, alpha.div, stat=med.stat, trace = TRUE, plot = FALSE)
#min.div<-zonal.stats(poly, alpha.div, stat=min.stat, trace = TRUE, plot = FALSE)

#mean.div50<-zonal.stats(poly, alpha.div50, stat=mean.stat, trace = TRUE, plot = FALSE)
#med.div50<-zonal.stats(poly, alpha.div50, stat=med.stat, trace = TRUE, plot = FALSE)
#min.div50<-zonal.stats(poly, alpha.div50, stat=min.stat, trace = TRUE, plot = FALSE)

#mean.div70<-zonal.stats(poly, alpha.div70, stat=mean.stat, trace = TRUE, plot = FALSE)
#med.div70<-zonal.stats(poly, alpha.div70, stat=med.stat, trace = TRUE, plot = FALSE)
#min.div70<-zonal.stats(poly, alpha.div70, stat=min.stat, trace = TRUE, plot = FALSE)


#length(mean.div)
#df$mean.div<-mean.div
#df$med.div<-med.div
#df$min.div<-min.div

#df$mean.div50<-mean.div50
#df$med.div50<-med.div50
#df$min.div50<-min.div50

#df$mean.div70<-mean.div70
#df$med.div70<-med.div70
#df$min.div70<-min.div70


