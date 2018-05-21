require(lubridate)
library(lmodel2)
library(raster)
library(rgdal)

source(file='/usr3/graduate/emelaas/Code/R/landsat_sentinel/v1_3/MCD12Q2C6_functions_1.R')

#List of study site locations
setwd('/projectnb/modislc/projects/landsat_sentinel/sites/')
d <- dir()
site_info <- rbind(c('COW','coweeta',35.0597,-83.4280,17,'17SKU'),
  c('HBB','hubbard',43.9438,-71.701,18,'18TYP'),
  c('HVF','harvard',42.5378,-72.1715,18,'18TYN'),
  c('UIE','uiefmaize',40.062822,-88.196128,16,'16TCK'))

#Generate table with headers
site_info <- data.frame(site_info[,1],site_info[,2],
  as.numeric(site_info[,3]),as.numeric(site_info[,4]),as.numeric(site_info[,5]),site_info[,6])
colnames(site_info) <- c('code','name','lat','lon','UTM','granule')

site <- site_info[3,]

pcam_proj <- readOGR('/projectnb/modislc/projects/landsat_sentinel/shapefiles/pcam_sites.shp')

deg2rad <- function(deg) {(deg * pi) / (180)}

slope <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site$code,'_slope_crop_wshed.tif',sep=''))
DEM <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site$code,'_DEM_crop_wshed.tif',sep=''))
conif <- raster('/projectnb/modislc/projects/landsat_sentinel/orthophoto/COW/conif_perc.tif')

aspect <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site$code,'_aspect_crop_wshed.tif',sep=''))
aspect_trans <- cos(deg2rad(aspect))

conif_vals <- getValues(conif)
aspect_vals <- getValues(aspect_trans)
aspect_vals_plot <- round(aspect_vals*10)/10
DEM_vals <- getValues(DEM)
slope_vals <- getValues(slope)

red_sd <- matrix(NA,240,11)
nir_sd <- matrix(NA,240,11)

all_meth <- c('none','cosine','improvedcosine','minnaert','minslope','ccorrection','ccorrection_smooth',
  'empirical','rotational','gamma','SCS','illumination')
#x11(h=8,w=11)
pdf(h=8,w=11,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/illumination.pdf')
par(mfrow=c(3,4),mar=c(2,2,2,2))
for (i in 1:11){
  print(i)
  
  load(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/',site$name,'_bands_topocorr_method',i,sep=''))
  
  time <- as.Date(b4[,4])
  yr <- b4[,1]
  doy <- b4[,2]
  sat <- b4[,3]
  index <- data.frame(yr,doy,sat,time)
  
  blue <- b1[,7:ncol(b1)]
  green <- b2[,7:ncol(b2)]
  red <- b3[,7:ncol(b3)]
  nir <- b4[,7:ncol(b4)]
  swir <- b5[,7:ncol(b5)]
  
  w1 <- which(aspect_vals >= 0.9 & conif_vals < 0.05) #North-facing pixels
  w2 <- which(aspect_vals <= -0.9 & conif_vals < 0.05) #South-facing pixels
  w3 <- which(DEM_vals <= 700) #Low elevation
  w4 <- which(DEM_vals > 1300) #High elevation
  w5 <- which(slope_vals<quantile(slope_vals,0.1) & conif_vals < 0.05)
  w6 <- which(slope_vals>quantile(slope_vals,0.9) & conif_vals < 0.05)
  wL <- which(sat=='L')
  wS <- which(sat=='S')
  wS_summer <- which(sat=='S' & yr==2016 & doy>180 & doy<260)
  wL_summer <- which(sat=='L' & yr==2016 & doy>180 & doy<260)
  wconif <- which(conif_vals<0.05)
  
  red_sd[,i] <- apply(red[,wconif],1,sd,na.rm=TRUE)
  nir_sd[,i] <- apply(nir[,wconif],1,sd,na.rm=TRUE)
  
  if (i == 11) plot(setValues(aspect,as.numeric(nir[34,])/10000),axes=FALSE,zlim=c(0,1),colNA='black',legend.width=1.2,
    main=all_meth[i])
  if (i != 11) plot(setValues(aspect,as.numeric(nir[34,])/10000),axes=FALSE,zlim=c(0,1),colNA='black',legend=FALSE,
    main=all_meth[i])
}
dev.off()

#### Topo Correction DOY vs Band comparison ####

#x11(h=9,w=11)
pdf(h=9,w=11,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/topocorr_comp_doy_nir.pdf')
par(mfrow=c(3,4),mar=c(4,4,3,2))
for (i in 1:11){
  plot(doy,nir_sd[,i]/10000-nir_sd[,1]/10000,ylim=c(-0.05,0.05),pch=16,cex=1.5,main=all_meth[i],
    xlab='DOY',ylab='Change in NIR StDev')
  abline(h=0,lty=2)
}
dev.off()

#x11(h=9,w=11)
pdf(h=9,w=11,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/topocorr_comp_doy_red.pdf')
par(mfrow=c(3,4),mar=c(4,4,3,2))
for (i in 1:11){
  plot(doy,red_sd[,i]/10000-red_sd[,1]/10000,ylim=c(-0.025,0.025),pch=16,cex=1.5,main=all_meth[i],
    xlab='DOY',ylab='Change in Red StDev')
  abline(h=0,lty=2)
}
dev.off()

#### Average Interval Greenup ####

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_gup_avg*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
gup_avg <- stack(in_dirs_tile)
gup_avg_vals <- getValues(gup_avg)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2_L8only/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_gup_avg*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
gup_avg2 <- stack(in_dirs_tile)
gup_avg2_vals <- getValues(gup_avg2)

my.colors = colorRampPalette(c("red","white","blue"))
#x11(h=8.2,w=9.4)
pdf(h=8.2,w=9.4,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/gup_avg_plot.pdf')
par(mfrow=c(2,3),mar=c(2,2,3,3))
for(i in 1:5){
  if (i==1)  plot(gup_avg2[[i]],axes=FALSE,zlim=c(0,50),col=my.colors(150),legend.width=1.2,main=(i+2012))
  if (i>=2 & i<=4)  plot(gup_avg2[[i]],axes=FALSE,zlim=c(0,50),col=my.colors(150),legend=FALSE,main=(i+2012))
  if (i==5)  plot(gup_avg[[4]],axes=FALSE,zlim=c(0,50),col=my.colors(150),legend=FALSE,main='2016 (HLS)')
}
boxplot(cbind(gup_avg_vals[,1:3],gup_avg2_vals[,4],gup_avg_vals[,4]),outline=FALSE,
  names=c('2013','2014','2015','2016','2016 (HLS)'),main='Greenup Average Obs. Interval (days)')
dev.off()


#### Longest Interval Greenup ####

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_gup_longest*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
gup_lon <- stack(in_dirs_tile)
gup_lon <- crop(gup_lon,pcam_proj)
gup_lon_vals <- getValues(gup_lon)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2_L8only/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_gup_longest*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
gup_lon2 <- stack(in_dirs_tile)
gup_lon2 <- crop(gup_lon2,pcam_proj)
gup_lon2_vals <- getValues(gup_lon2)

my.colors = colorRampPalette(c("red","white","blue"))
x11(h=6.2,w=9.4)
#pdf(h=8.2,w=9.4,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/gup_lon_plot.pdf')
par(mfrow=c(2,3),mar=c(2,2,3,3))
for(i in 1:5){
  if (i==1)  plot(gup_lon2[[i]],axes=FALSE,zlim=c(0,120),col=my.colors(150),legend.width=1.2,main=(i+2012))
  if (i>=2 & i<=4)  plot(gup_lon2[[i]],axes=FALSE,zlim=c(0,120),col=my.colors(150),legend=FALSE,main=(i+2012))
  if (i==5)  plot(gup_lon[[4]],axes=FALSE,zlim=c(0,120),col=my.colors(150),legend=FALSE,main='2016 (HLS)')
  points(pcam_proj,pch=21,col='black',bg='yellow',cex=1.5)
}
boxplot(cbind(gup_lon_vals[,1:3],gup_lon2_vals[,4],gup_lon_vals[,4]),outline=FALSE,
  names=c('2013','2014','2015','2016','2016 (HLS)'),main='Greenup Longest Obs. Interval (days)')
dev.off()


#### Average Interval Greendown ####

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_gdown_avg*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
gdown_avg <- stack(in_dirs_tile)
gdown_avg_vals <- getValues(gdown_avg)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2_L8only/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_gdown_avg*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
gdown_avg2 <- stack(in_dirs_tile)
gdown_avg2_vals <- getValues(gdown_avg2)

my.colors = colorRampPalette(c("red","white","blue"))
#x11(h=6.2,w=9.4)
pdf(h=8.2,w=9.4,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/gdown_avg_plot.pdf')
par(mfrow=c(2,3),mar=c(2,2,3,3))
plot(gdown_avg2[[2]],axes=FALSE,zlim=c(0,30),col=my.colors(150),legend.width=1.2,main='2014')
plot(gdown_avg2[[3]],axes=FALSE,zlim=c(0,30),col=my.colors(150),legend=FALSE,main='2015')
plot(gdown_avg[[3]],axes=FALSE,zlim=c(0,30),col=my.colors(150),legend=FALSE,main='2015 (HLS)')
plot(gdown_avg2[[4]],axes=FALSE,zlim=c(0,30),col=my.colors(150),legend=FALSE,main='2016')
plot(gdown_avg[[4]],axes=FALSE,zlim=c(0,30),col=my.colors(150),legend=FALSE,main='2016 (HLS)')
boxplot(cbind(gdown_avg2_vals[,2:3],gdown_avg_vals[,3],gdown_avg2_vals[,4],gdown_avg_vals[,4]),outline=FALSE,
  names=c('2014','2015','2015 (HLS)','2016','2016 (HLS)'),main='Senescence Average Obs. Interval (days)')
dev.off()


#### Longest Interval Greendown ####

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_gdown_lon*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
gdown_lon <- stack(in_dirs_tile)
gdown_lon_vals <- getValues(gdown_lon)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2_L8only/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_gdown_lon*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
gdown_lon2 <- stack(in_dirs_tile)
gdown_lon2_vals <- getValues(gdown_lon2)

my.colors = colorRampPalette(c("red","white","blue"))
#x11(h=6.2,w=9.4)
pdf(h=8.2,w=9.4,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/gdown_lon_plot.pdf')
par(mfrow=c(2,3),mar=c(2,2,3,3))
plot(gdown_lon2[[2]],axes=FALSE,zlim=c(0,100),col=my.colors(150),legend.width=1.2,main='2014')
plot(gdown_lon2[[3]],axes=FALSE,zlim=c(0,100),col=my.colors(150),legend=FALSE,main='2015')
plot(gdown_lon[[3]],axes=FALSE,zlim=c(0,100),col=my.colors(150),legend=FALSE,main='2015 (HLS)')
plot(gdown_lon2[[4]],axes=FALSE,zlim=c(0,100),col=my.colors(150),legend=FALSE,main='2016')
plot(gdown_lon[[4]],axes=FALSE,zlim=c(0,100),col=my.colors(150),legend=FALSE,main='2016 (HLS)')
boxplot(cbind(gdown_lon2_vals[,2:3],gdown_lon_vals[,3],gdown_lon2_vals[,4],gdown_lon_vals[,4]),outline=FALSE,
  names=c('2014','2015','2015 (HLS)','2016','2016 (HLS)'),main='Senescence Longest Obs. Interval (days)')
dev.off()


#### Start of Spring ####

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_SOS*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
SOS <- stack(in_dirs_tile)
SOS_vals <- getValues(SOS)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_1_SOS*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
SOS1 <- stack(in_dirs_tile)
SOS1_vals <- getValues(SOS1)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2_L8only/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_SOS*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
SOS2 <- stack(in_dirs_tile)
SOS2_vals <- getValues(SOS2)

my.colors = colorRampPalette(c("red","yellow","blue"))
#x11(h=6.2,w=9.4)
pdf(h=8.2,w=9.4,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/SOS_plot.pdf')
par(mfrow=c(2,3),mar=c(2,2,3,3))
plot(SOS2[[1]],axes=FALSE,zlim=c(30,130),col=my.colors(150),legend.width=1.2,main='2013',colNA='black')
plot(SOS2[[2]],axes=FALSE,zlim=c(30,130),col=my.colors(150),legend=FALSE,main='2014',colNA='black')
plot(SOS2[[3]],axes=FALSE,zlim=c(30,130),col=my.colors(150),legend=FALSE,main='2015',colNA='black')
plot(SOS2[[4]],axes=FALSE,zlim=c(30,130),col=my.colors(150),legend=FALSE,main='2016',colNA='black')
plot(SOS[[4]],axes=FALSE,zlim=c(30,130),col=my.colors(150),legend=FALSE,main='2016 (HLS)',colNA='black')
boxplot(cbind(SOS2_vals[,1:4],SOS_vals[,4]),outline=FALSE,
  names=c('2013','2014','2015','2016','2016 (HLS)'),main='Start of Spring (DOY)')
dev.off()

my.colors = colorRampPalette(c("red","yellow","blue"))
#x11(h=5,w=11)
pdf(h=5,w=11,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/SOS_plot_2016.pdf')
par(mfrow=c(1,4),mar=c(2,2,3,3))
plot(SOS1[[4]],axes=FALSE,zlim=c(30,130),col=my.colors(150),legend.width=1.2,main='No corr. (HLS)',colNA='black')
plot(SOS2[[4]],axes=FALSE,zlim=c(30,130),col=my.colors(150),legend=FALSE,main='Corr. (L8 only)',colNA='black')
plot(SOS[[4]],axes=FALSE,zlim=c(30,130),col=my.colors(150),legend=FALSE,main='Corr. (HLS)',colNA='black')
boxplot(cbind(SOS1_vals[,4],SOS2_vals[,4],SOS_vals[,4]),outline=FALSE,
  names=c('Uncorr. (HLS)','Corr. (L8)','Corr. (HLS)'),main='2016 Start of Spring (DOY)')
dev.off()

#### Middle of Spring ####

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_MOS*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
MOS <- stack(in_dirs_tile)
MOS <- crop(MOS,pcam_proj)
MOS_vals <- getValues(MOS)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_1_MOS*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
MOS1 <- stack(in_dirs_tile)
MOS1 <- crop(MOS1,pcam_proj)
MOS1_vals <- getValues(MOS1)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2_L8only/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_MOS*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
MOS2 <- stack(in_dirs_tile)
MOS2 <- crop(MOS2,pcam_proj)
MOS2_vals <- getValues(MOS2)

my.colors = colorRampPalette(c("red","yellow","blue"))
#x11(h=6.2,w=9.4)
pdf(h=8.2,w=9.4,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/MOS_plot.pdf')
par(mfrow=c(2,3),mar=c(2,2,3,3))
plot(MOS2[[1]],axes=FALSE,zlim=c(90,160),col=my.colors(150),legend.width=1.2,main='2013',colNA='black')
plot(MOS2[[2]],axes=FALSE,zlim=c(90,160),col=my.colors(150),legend=FALSE,main='2014',colNA='black')
plot(MOS2[[3]],axes=FALSE,zlim=c(90,160),col=my.colors(150),legend=FALSE,main='2015',colNA='black')
plot(MOS2[[4]],axes=FALSE,zlim=c(90,160),col=my.colors(150),legend=FALSE,main='2016',colNA='black')
plot(MOS[[4]],axes=FALSE,zlim=c(90,160),col=my.colors(150),legend=FALSE,main='2016 (HLS)',colNA='black')
boxplot(cbind(MOS2_vals[,1:4],MOS_vals[,4]),outline=FALSE,
  names=c('2013','2014','2015','2016','2016 (HLS)'),main='Middle of Spring (DOY)')
dev.off()

my.colors = colorRampPalette(c("red","orange","yellow","green","blue"))
x11(h=4.75,w=13)
#pdf(h=4.75,w=13,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/HVF/MOS_plot_2016.pdf')
par(mfrow=c(1,4),mar=c(2,2,3,3))
plot(MOS1[[4]],axes=FALSE,zlim=c(105,150),col=my.colors(150),legend.width=1.4,main='No corr. (HLS)',colNA='black')
points(pcam_proj,cex=1.75,pch=21,col='black',bg='yellow')
plot(MOS2[[4]],axes=FALSE,zlim=c(105,150),col=my.colors(150),legend=FALSE,main='Corr. (L8 only)',colNA='black')
points(pcam_proj,cex=1.75,pch=21,col='black',bg='yellow')
plot(MOS[[4]],axes=FALSE,zlim=c(105,150),col=my.colors(150),legend=FALSE,main='Corr. (HLS)',colNA='black')
points(pcam_proj,cex=1.75,pch=21,col='black',bg='yellow')
boxplot(cbind(MOS1_vals[,4],MOS2_vals[,4],MOS_vals[,4]),outline=FALSE,
  names=c('Uncorr. (HLS)','Corr. (L8)','Corr. (HLS)'),main='2016 Middle of Spring (DOY)')
dev.off()

my.colors = colorRampPalette(c("red","orange","yellow","green","blue"))
#x11(h=4.75,w=13)
pdf(h=4.75,w=13,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/HVF/MOS_plot_2016.pdf')
par(mfrow=c(1,4),mar=c(2,2,3,3))
plot(MOS1[[4]],axes=FALSE,zlim=c(105,150),main='No correction (HLS)',legend=FALSE,col=my.colors(150),colNA='black')
plot(MOS1[[4]],zlim=c(105,150),legend.only=T,legend.shrink=1,col=my.colors(150),legend.width=1.2,
  axis.args=list(at=seq(105,150,10),labels=seq(105,150,10),cex.axis=1.3,font=2))
points(pcam_proj,cex=2.75,pch=4,col='black',bg='yellow')
plot(MOS2[[4]],axes=FALSE,zlim=c(105,150),col=my.colors(150),legend=FALSE,main='Corrected (L8)',cex.main=1.3,colNA='black')
points(pcam_proj,cex=2.75,pch=4,col='black',bg='yellow')
plot(MOS[[4]],axes=FALSE,zlim=c(105,150),col=my.colors(150),legend=FALSE,main='Corrected (HLS)',cex.main=1.3,colNA='black')
points(pcam_proj,cex=2.75,pch=4,col='black',bg='yellow')
par(font.axis=2)
boxplot(cbind(MOS1_vals[,4],MOS2_vals[,4],MOS_vals[,4]),outline=FALSE,cex.main=1.3,
  names=c('No (HLS)','Corr (L8)','Corr (HLS)'),main='2016 Middle of Spring (DOY)',cex.axis=1.3)
dev.off()


#### End of Spring ####

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_EOS*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
EOS <- stack(in_dirs_tile)
EOS_vals <- getValues(EOS)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_1_EOS*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
EOS1 <- stack(in_dirs_tile)
EOS1_vals <- getValues(EOS1)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2_L8only/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_EOS*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
EOS2 <- stack(in_dirs_tile)
EOS2_vals <- getValues(EOS2)

my.colors = colorRampPalette(c("red","yellow","blue"))
#x11(h=6.2,w=9.4)
pdf(h=8.2,w=9.4,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/EOS_plot.pdf')
par(mfrow=c(2,3),mar=c(2,2,3,3))
plot(EOS2[[1]],axes=FALSE,zlim=c(130,200),col=my.colors(150),legend.width=1.2,main='2013',colNA='black')
plot(EOS2[[2]],axes=FALSE,zlim=c(130,200),col=my.colors(150),legend=FALSE,main='2014',colNA='black')
plot(EOS2[[3]],axes=FALSE,zlim=c(130,200),col=my.colors(150),legend=FALSE,main='2015',colNA='black')
plot(EOS2[[4]],axes=FALSE,zlim=c(130,200),col=my.colors(150),legend=FALSE,main='2016',colNA='black')
plot(EOS[[4]],axes=FALSE,zlim=c(130,200),col=my.colors(150),legend=FALSE,main='2016 (HLS)',colNA='black')
boxplot(cbind(EOS2_vals[,1:4],EOS_vals[,4]),outline=FALSE,
  names=c('2013','2014','2015','2016','2016 (HLS)'),main='End of Spring (DOY)')
dev.off()

my.colors = colorRampPalette(c("red","yellow","blue"))
#x11(h=6.8,w=7.2)
pdf(h=5,w=11,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/EOS_plot_2016.pdf')
par(mfrow=c(1,4),mar=c(2,2,3,3))
plot(EOS1[[4]],axes=FALSE,zlim=c(130,200),col=my.colors(150),legend.width=1.2,main='No corr. (HLS)',colNA='black')
plot(EOS2[[4]],axes=FALSE,zlim=c(130,200),col=my.colors(150),legend=FALSE,main='Corr. (L8 only)',colNA='black')
plot(EOS[[4]],axes=FALSE,zlim=c(130,200),col=my.colors(150),legend=FALSE,main='Corr. (HLS)',colNA='black')
boxplot(cbind(EOS1_vals[,4],EOS2_vals[,4],EOS_vals[,4]),outline=FALSE,
  names=c('Uncorr. (HLS)','Corr. (L8)','Corr. (HLS)'),main='2016 End of Spring (DOY)')
dev.off()

#### Middle of Autumn ####

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_MOA*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
MOA <- stack(in_dirs_tile)
MOA_vals <- getValues(MOA)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_1_MOA*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
MOA1 <- stack(in_dirs_tile)
MOA1_vals <- getValues(MOA1)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2_L8only/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_MOA*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
MOA2 <- stack(in_dirs_tile)
MOA2_vals <- getValues(MOA2)

my.colors = colorRampPalette(c("red","yellow","blue"))
#x11(h=6.2,w=9.4)
pdf(h=8.2,w=9.4,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/MOA_plot.pdf')
par(mfrow=c(2,3),mar=c(2,2,3,3))
plot(MOA2[[2]],axes=FALSE,zlim=c(260,330),col=my.colors(150),legend.width=1.2,main='2014',colNA='black')
plot(MOA2[[3]],axes=FALSE,zlim=c(260,330),col=my.colors(150),legend=FALSE,main='2015',colNA='black')
plot(MOA[[3]],axes=FALSE,zlim=c(260,330),col=my.colors(150),legend=FALSE,main='2015 (HLS)',colNA='black')
plot(MOA2[[4]],axes=FALSE,zlim=c(260,330),col=my.colors(150),legend=FALSE,main='2016',colNA='black')
plot(MOA[[4]],axes=FALSE,zlim=c(260,330),col=my.colors(150),legend=FALSE,main='2016 (HLS)',colNA='black')
boxplot(cbind(MOA2_vals[,2:3],MOA_vals[,3],MOA2_vals[,4],MOA_vals[,4]),outline=FALSE,
  names=c('2014','2015','2015 (HLS)','2016','2016 (HLS)'),main='Middle of Autumn (DOY)')
dev.off()

my.colors = colorRampPalette(c("red","yellow","blue"))
#x11(h=6.8,w=7.2)
pdf(h=5,w=11,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/MOA_plot_2016.pdf')
par(mfrow=c(1,4),mar=c(2,2,3,3))
plot(MOA1[[4]],axes=FALSE,zlim=c(260,330),col=my.colors(150),legend.width=1.2,main='No corr. (HLS)',colNA='black')
plot(MOA2[[4]],axes=FALSE,zlim=c(260,330),col=my.colors(150),legend=FALSE,main='Corr. (L8 only)',colNA='black')
plot(MOA[[4]],axes=FALSE,zlim=c(260,330),col=my.colors(150),legend=FALSE,main='Corr. (HLS)',colNA='black')
boxplot(cbind(MOA1_vals[,4],MOA2_vals[,4],MOA_vals[,4]),outline=FALSE,
  names=c('Uncorr. (HLS)','Corr. (L8)','Corr. (HLS)'),main='2016 End of Spring (DOY)')
dev.off()


#### Single Pixel Plot ####
source(file='/usr3/graduate/emelaas/Code/R/landsat_sentinel/v1_3/MCD12Q2C6_functions_1.R')
load(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/',site$name,'_bands_topocorr_method',9,sep=''))

time <- as.Date(b4[,4])
yr <- b4[,1]
doy <- b4[,2]
sat <- b4[,3]
index <- data.frame(yr,doy,sat,time)

ts1 <- evi2[,20464]/10000
psl1 <- PlotSeriesLandsat3(ts1[which(yr>=2015 & yr<=2016)],time[which(yr>=2015 & yr<=2016)],pheno_pars,
  'Golf Course','EVI2')
abline(v=psl1[[3]][[1]][[2]],lty=3,col='darkgreen')
abline(v=psl1[[3]][[1]][[5]],lty=3,col='darkred')
points(time[which(sat=='S')],ts1[which(sat=='S')],pch=21,col='black',bg='yellow',cex=1)
abline(v=tmp[[3]][[1]][[2]][1:2],lwd=2)
abline(v=tmp[[3]][[1]][[5]][1:2],lwd=2)

ts2 <- evi2[,11604]/10000
psl2 <- PlotSeriesLandsat3(ts2[which(yr>=2015 & yr<=2016)],time[which(yr>=2015 & yr<=2016)],pheno_pars,
  'Deciduous Broadleaf Forest','EVI2')
abline(v=psl2[[3]][[1]][[2]],lty=3,col='darkgreen')
abline(v=psl2[[3]][[1]][[5]],lty=3,col='darkred')
points(time[which(sat=='S')],ts2[which(sat=='S')],pch=21,col='black',bg='yellow',cex=1)
abline(v=tmp[[3]][[1]][[2]][2],lwd=2)
abline(v=tmp[[3]][[1]][[5]][2],lwd=2)

ts2 <- evi2[,9371]/10000
psl2 <- PlotSeriesLandsat3(ts2[which(yr>=2014 & yr<=2016)],time[which(yr>=2014 & yr<=2016)],pheno_pars,
  'Evergreen Needleleaf Forest','EVI2')
ts3 <- evi2[,20650]/10000
psl2 <- PlotSeriesLandsat3(ts3[which(yr>=2014 & yr<=2016)],time[which(yr>=2014 & yr<=2016)],pheno_pars,
  'Grassland','EVI2')

ts1 <- evi2[,7503]/10000
ts2 <- evi2[,9371]/10000
ts3 <- evi2[,20650]/10000
df <- data.frame(yr,doy,sat,ts1,ts2,ts3)
#x11()
pdf('/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/HVF/pix_sample_ts.pdf')
par(mfrow=c(3,1),mar=c(3,4,3,3),font.axis=2,font.lab=2)
plot(time,ts1,pch=16,cex=1.3,main='Deciduous Broadleaf Forest',ylab='EVI2',cex.axis=1.2,cex.main=1.3,cex.lab=1.2)
points(time[wS],ts1[wS],pch=21,col='black',bg='yellow',cex=1.3)
legend('topleft',c('S2A','L8'),pch=c(21,16),pt.bg=c('yellow','black'),bty='n',cex=1.2,text.font=2)
plot(time,ts2,pch=16,cex=1.3,main='Evergreen Needleleaf Forest',ylab='EVI2',cex.axis=1.2,cex.main=1.3,cex.lab=1.2)
points(time[wS],ts2[wS],pch=21,col='black',bg='yellow',cex=1.3)
plot(time,ts3,pch=16,cex=1.3,main='Grassland',ylab='EVI2',cex.axis=1.2,cex.main=1.3,cex.lab=1.2)
points(time[wS],ts3[wS],pch=21,col='black',bg='yellow',cex=1.3)
dev.off()

pix <- 7839
ts1 <- red[,pix]/10000
ts2 <- nir[,pix]/10000
ts3 <- ndvi[,pix]/10000
ts4 <- evi2[,pix]/10000
PlotSeriesLandsat2(ts4,time,pheno_pars,'coweeta-8672','EVI2')



#### In situ observations ####
load('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/harvard_pcam_dates')

#HVF
sos50_sat <- t(extract(MOS[[2:4]],pcam_proj))
eos50_sat <- t(extract(MOA[[2:4]],pcam_proj))

#x11(h=5,w=9)
pdf(h=5,w=9,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/HVF/insitu_vs_HLS.pdf')
par(mfrow=c(1,2))
plot(as.numeric(sos50),as.numeric(sos50_sat),xlim=c(100,150),ylim=c(100,150),pch=16,col='blue',
  ylab='HLS MOS (DOY)',xlab='In Situ MOS (DOY)',cex=1.5)
abline(0,1,lty=2)
legend('topleft',legend='PhenoCam',pch=16,col='blue',bty='n')

y1 <- as.numeric(sos50_sat)
x1 <- as.numeric(sos50)
lm2_1 <- lmodel2(y1~x1)
r2_1 <- lm2_1$rsquare
cortest <- cor.test(x1,y1,use='pairwise.complete.obs',method='pearson')
my.p_1 <- cortest$p.value
rmse_1 <- sqrt(mean((y1-x1)^2,na.rm=TRUE))
rp = vector('expression',2)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) < 0.01,')')),
  list(MYVALUE1 = format(r2_1,dig=2),MYVALUE2 = format(my.p_1,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse_1,dig=2)))[2]
legend(125,110, legend = rp, bty = 'n',cex=0.9)

plot(as.numeric(eos50),as.numeric(eos50_sat),xlim=c(100,365),ylim=c(100,365),pch=16,col='blue',
  ylab='HLS MOA (DOY)',xlab='In Situ MOA (DOY)',cex=1.5)
abline(0,1,lty=2)
y1 <- as.numeric(eos50_sat)
x1 <- as.numeric(eos50)
lm2_1 <- lmodel2(y1~x1)
r2_1 <- lm2_1$rsquare
cortest <- cor.test(x1,y1,use='pairwise.complete.obs',method='pearson')
my.p_1 <- cortest$p.value
rmse_1 <- sqrt(mean((y1-x1)^2,na.rm=TRUE))
rp = vector('expression',2)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) == MYVALUE2,')')),
  list(MYVALUE1 = format(r2_1,dig=2),MYVALUE2 = format(my.p_1,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse_1,dig=2)))[2]
legend(230,150, legend = rp, bty = 'n',cex=0.9)

dev.off()

