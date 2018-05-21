library(lmodel2)
library(raster)
library(rgdal) 

source(file='/usr3/graduate/emelaas/Code/R/landsat_sentinel/v1_3/MCD12Q2C6_functions_1.R')

#List of study site locations
setwd('/projectnb/modislc/projects/landsat_sentinel/sites/')
d <- dir()
site_info <- rbind(c('COW','coweeta',35.0597,-83.4280,17,'17SKU'),
  c('HBB','hubbard',43.9438,-71.701,18,'18TYP'))

#Generate table with headers
site_info <- data.frame(site_info[,1],site_info[,2],
  as.numeric(site_info[,3]),as.numeric(site_info[,4]),as.numeric(site_info[,5]),site_info[,6])
colnames(site_info) <- c('code','name','lat','lon','UTM','granule')

site <- site_info[1,]

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

my.colors1 = colorRampPalette(c("white","yellow","orange","black"))
my.colors2 = colorRampPalette(c("gray25","gray50","gray75","white"))

all_meth <- c('None','Cosine','Improved Cosine','Minnaert','Minnaert w/ Slope','C Correction','Smoothed C Correction',
  'Empirical','Rotational','Gamma','SCS','Illumination')
x11(h=8,w=11)
#pdf(h=8,w=11,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/COW/topocorr_comp_2014014.pdf')
par(mfrow=c(3,4),mar=c(2,2,2,2))
plot(aspect_trans,col=my.colors1(150),axes=FALSE,legend.width=1.4,main='aspect')
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
  
  if (i == 1) plot(setValues(aspect,as.numeric(nir[34,])/10000),axes=FALSE,zlim=c(0,0.5),colNA='black',legend.width=1.3,
    main=all_meth[i],col=my.colors2(150))
  if (i != 1) plot(setValues(aspect,as.numeric(nir[34,])/10000),axes=FALSE,zlim=c(0,0.5),colNA='black',legend=FALSE,
    main=all_meth[i],col=my.colors2(150))
}
dev.off()


my.colors1 = colorRampPalette(c("white","yellow","orange","black"))
my.colors2 = colorRampPalette(c("gray25","gray50","gray75","white"))

#x11(h=4.2,w=13)
pdf(h=4.2,w=13,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/COW/topocorr_comp2_2014014.pdf')
par(mfrow=c(1,3),mar=c(2,2,2,2))
plot(aspect_trans,col=my.colors1(150),axes=FALSE,legend=FALSE,main='Transformed Aspect',cex.main=1.2)
plot(aspect_trans,zlim=c(-1,1),legend.only=T,legend.shrink=1,col=my.colors1(150),legend.width=1.75,
  axis.args=list(at=seq(-1,1,0.5),labels=seq(-1,1,0.5),cex.axis=1.4,font=2))
for (i in c(1,9)){
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
    
  if (i == 1) {
    plot(setValues(aspect,as.numeric(nir[34,])/10000),axes=FALSE,zlim=c(0,0.5),colNA='black',legend=FALSE,
      main='Uncorrected NIR',col=my.colors2(150),cex.main=1.2)
    plot(aspect_trans,zlim=c(0,0.5),legend.only=T,legend.shrink=1,col=my.colors2(150),legend.width=1.75,
      axis.args=list(at=seq(0,0.5,0.1),labels=seq(0,0.5,0.1),cex.axis=1.4,font=2))
  } 
  if (i == 9) {
    plot(setValues(aspect,as.numeric(nir[34,])/10000),axes=FALSE,zlim=c(0,0.5),colNA='black',legend=FALSE,
      main='Rotational Correction',col=my.colors2(150),cex.main=1.2)
#     plot(aspect_trans,zlim=c(0,0.5),legend.only=T,legend.shrink=1,col=my.colors2(150),legend.width=1.75,
#       axis.args=list(at=seq(0,0.5,0.1),labels=seq(0,0.5,0.1),cex.axis=1.3,font=2))
  } 
}
dev.off()

#### Topo Correction DOY vs Band comparison ####

#x11(h=9,w=11)
pdf(h=9,w=11,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/COW/topocorr_comp_doy_nir.pdf')
par(mfrow=c(3,4),mar=c(4,4,3,2))
for (i in 1:11){
  plot(doy,nir_sd[,i]/10000-nir_sd[,1]/10000,ylim=c(-0.05,0.05),pch=16,cex=1.5,main=all_meth[i],
    xlab='DOY',ylab='Change in NIR StDev')
  abline(h=0,lty=2)
}
dev.off()

#x11(h=9,w=11)
pdf(h=9,w=11,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/COW/topocorr_comp_doy_red.pdf')
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
#x11(h=3.45,w=14)
pdf(h=3.45,w=14,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/COW/gup_avg_plot.pdf')
par(mfrow=c(1,4),mar=c(2,2,3,3))
for(i in 3:5){
  if (i==3) {
    r <- gup_avg2[[i]]
    plot(r,axes=FALSE,zlim=c(0,50),main='2015',legend=FALSE,col=my.colors(150))
    plot(r,zlim=c(0,50),legend.only=T,legend.shrink=1,col=my.colors(150),legend.width=1.2,
      axis.args=list(at=seq(0,50,10),labels=seq(0,50,10),cex.axis=1.3,font=2))
  }  
  if (i==4)  plot(gup_avg2[[i]],axes=FALSE,zlim=c(0,50),col=my.colors(150),legend=FALSE,main='2016 (L8)',cex.main=1.3)
  if (i==5)  plot(gup_avg[[4]],axes=FALSE,zlim=c(0,50),col=my.colors(150),legend=FALSE,main='2016 (HLS)',cex.main=1.3)
}
par(font.axis=2)
b <- boxplot(cbind(gup_avg_vals[,3],gup_avg2_vals[,4],gup_avg_vals[,4]),outline=FALSE,cex.main=1.3,
  names=c('2015','2016 (L8)','2016 (HLS)'),main='Greenup Average Obs. Interval (days)',cex.axis=1.3)
dev.off()


#### Longest Interval Greenup ####

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_gup_longest*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
gup_lon <- stack(in_dirs_tile)
gup_lon_vals <- getValues(gup_lon)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2_L8only/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_gup_longest*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
gup_lon2 <- stack(in_dirs_tile)
gup_lon2_vals <- getValues(gup_lon2)

my.colors = colorRampPalette(c("red","white","blue"))
#x11(h=3.45,w=14)
pdf(h=3.45,w=14,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/COW/gup_lon_plot.pdf')
par(mfrow=c(1,4),mar=c(2,2,3,3))
for(i in 3:5){
  if (i==3) {
    r <- gup_lon2[[i]]
    plot(r,axes=FALSE,zlim=c(0,120),main='2015',legend=FALSE,col=my.colors(150))
    plot(r,zlim=c(0,120),legend.only=T,legend.shrink=1,col=my.colors(150),legend.width=1.2,
      axis.args=list(at=seq(0,120,20),labels=seq(0,120,20),cex.axis=1.3,font=2))
  }  
  if (i==4)  plot(gup_lon2[[i]],axes=FALSE,zlim=c(0,120),col=my.colors(150),legend=FALSE,main='2016 (L8)',cex.main=1.3)
  if (i==5)  plot(gup_lon[[4]],axes=FALSE,zlim=c(0,120),col=my.colors(150),legend=FALSE,main='2016 (HLS)',cex.main=1.3)
}
par(font.axis=2)
boxplot(cbind(gup_lon_vals[,3],gup_lon2_vals[,4],gup_lon_vals[,4]),outline=FALSE,cex.main=1.3,
  names=c('2015','2016 (L8)','2016 (HLS)'),main='Greenup Longest Obs. Interval (days)',cex.axis=1.3)
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

my.colors = colorRampPalette(c("red","white","blue"))
#x11(h=3.45,w=14)
pdf(h=3.45,w=14,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/COW/gdown_avg_plot.pdf')
par(mfrow=c(1,4),mar=c(2,2,3,3))
for(i in 3:5){
  if (i==3) {
    r <- gdown_avg2[[i]]
    plot(r,axes=FALSE,zlim=c(0,30),main='2015',legend=FALSE,col=my.colors(150))
    plot(r,zlim=c(0,30),legend.only=T,legend.shrink=1,col=my.colors(150),legend.width=1.2,
      axis.args=list(at=seq(0,30,5),labels=seq(0,30,5),cex.axis=1.3,font=2))
  }  
  if (i==4)  plot(gdown_avg2[[i]],axes=FALSE,zlim=c(0,30),col=my.colors(150),legend=FALSE,main='2016 (L8)',cex.main=1.3)
  if (i==5)  plot(gdown_avg[[4]],axes=FALSE,zlim=c(0,30),col=my.colors(150),legend=FALSE,main='2016 (HLS)',cex.main=1.3)
}
par(font.axis=2)
b <- boxplot(cbind(gdown_avg_vals[,3],gdown_avg2_vals[,4],gdown_avg_vals[,4]),outline=FALSE,cex.main=1.3,
  names=c('2015','2016 (L8)','2016 (HLS)'),main='Greendown Average Obs. Interval (days)',cex.axis=1.3)
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

my.colors = colorRampPalette(c("blue","green","yellow","orange","red"))
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

my.colors = colorRampPalette(c("blue","green","yellow","orange","red"))
x11(h=2.75,w=11)
#pdf(h=5,w=11,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/SOS_plot_2016.pdf')
par(mfrow=c(1,4),mar=c(2,2,3,3))
plot(SOS1[[4]],axes=FALSE,zlim=c(30,130),col=my.colors(150),legend.width=1.4,main='No corr. (HLS)',colNA='black')
plot(SOS2[[4]],axes=FALSE,zlim=c(30,130),col=my.colors(150),legend=FALSE,main='Corr. (L8 only)',colNA='black')
plot(SOS[[4]],axes=FALSE,zlim=c(30,130),col=my.colors(150),legend=FALSE,main='Corr. (HLS)',colNA='black')
boxplot(cbind(SOS1_vals[,4],SOS2_vals[,4],SOS_vals[,4]),outline=FALSE,
  names=c('Uncorr. (HLS)','Corr. (L8)','Corr. (HLS)'),main='2016 Start of Spring (DOY)')
dev.off()

#### Middle of Spring ####

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_MOS*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
MOS <- stack(in_dirs_tile)
MOS_vals <- getValues(MOS)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_1_MOS*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
MOS1 <- stack(in_dirs_tile)
MOS1_vals <- getValues(MOS1)

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2_L8only/')
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste(site$code,'_9_MOS*.tif',sep='')),full.names=T,include.dirs=T,recursive=TRUE)
MOS2 <- stack(in_dirs_tile)
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

my.colors = colorRampPalette(c("red","orange","yellow","darkgreen","blue"))
#x11(h=3.45,w=14)
pdf(h=3.45,w=14,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/COW/MOS_plot_2016.pdf')
par(mfrow=c(1,4),mar=c(2,2,3,3))
plot(MOS1[[4]],axes=FALSE,zlim=c(90,160),main='No correction (HLS)',legend=FALSE,col=my.colors(150),colNA='black')
plot(MOS1[[4]],zlim=c(90,160),legend.only=T,legend.shrink=1,col=my.colors(150),legend.width=1.75,
  axis.args=list(at=seq(90,160,10),labels=seq(90,160,10),cex.axis=1.3,font=2))
plot(MOS2[[4]],axes=FALSE,zlim=c(90,160),col=my.colors(150),legend=FALSE,main='Corrected (L8)',cex.main=1.3,colNA='black')
plot(MOS[[4]],axes=FALSE,zlim=c(90,160),col=my.colors(150),legend=FALSE,main='Corrected (HLS)',cex.main=1.3,colNA='black')
par(font.axis=2)
b <- boxplot(cbind(MOS1_vals[,4],MOS2_vals[,4],MOS_vals[,4]),outline=FALSE,cex.main=1.3,
  names=c('No (HLS)','Corr (L8)','Corr (HLS)'),main='',cex.axis=1.3)
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
b <- boxplot(cbind(MOA2_vals[,2:3],MOA_vals[,3],MOA2_vals[,4],MOA_vals[,4]),outline=FALSE,
  names=c('2014','2015','2015 (HLS)','2016','2016 (HLS)'),main='Middle of Autumn (DOY)')
dev.off()

#my.colors = colorRampPalette(c("blue","darkgreen","yellow","orange","red"))
my.colors = colorRampPalette(c("red","orange","yellow","darkgreen","blue"))
#x11(h=3.45,w=14)
pdf(h=3.45,w=14,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/COW/MOA_plot_2016.pdf')
par(mfrow=c(1,4),mar=c(2,2,3,3))
plot(MOA1[[4]],axes=FALSE,zlim=c(260,330),main='No correction (HLS)',legend=FALSE,col=my.colors(150),colNA='black')
plot(MOA1[[4]],zlim=c(260,330),legend.only=T,legend.shrink=1,col=my.colors(150),legend.width=1.75,
  axis.args=list(at=seq(260,330,10),labels=seq(260,330,10),cex.axis=1.3,font=2))
plot(MOA2[[4]],axes=FALSE,zlim=c(260,330),col=my.colors(150),legend=FALSE,main='Corrected (L8)',cex.main=1.3,colNA='black')
plot(MOA[[4]],axes=FALSE,zlim=c(260,330),col=my.colors(150),legend=FALSE,main='Corrected (HLS)',cex.main=1.3,colNA='black')
par(font.axis=2)
b <- boxplot(cbind(MOA1_vals[,4],MOA2_vals[,4],MOA_vals[,4]),outline=FALSE,cex.main=1.3,
  names=c('No (HLS)','Corr (L8)','Corr (HLS)'),main='',cex.axis=1.3)
dev.off()


#### MOS vs. Elevation ####

#x11(h=7.2,w=7.2)
pdf(h=7.2,w=7.2,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/COW/MOS_MOA_DEM_aspect.pdf')
par(mfrow=c(2,2),mar=c(4,4,2,2))

MOS.lm.2016 <- lm(MOS1_vals[,4]~DEM_vals); s1 <- MOS.lm.2016$coefficients[2]
MOS.lm.2016a <- lm(MOS2_vals[,4]~DEM_vals); s2 <- MOS.lm.2016a$coefficients[2]
MOS.lm.2016b <- lm(MOS_vals[,4]~DEM_vals); s3 <- MOS.lm.2016b$coefficients[2]
#pdf(file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/MOS_vs_DEM.pdf')
plot(NA,xlim=c(600,1600),ylim=c(90,170),xlab='Elevation (m)',ylab='MOS (DOY)')
abline(MOS.lm.2016,col='blue',lwd=3)
abline(MOS.lm.2016a,col='black',lwd=3)
abline(MOS.lm.2016b,col='gray',lwd=3)
legend('topleft',
  c(paste('2016 (No Corr.) (',round(s1,3),')',sep=''),
    paste('2016 (L8 Only) (',round(s2,3),')',sep=''),
    paste('2016 (HLS) (',round(s3,3),')',sep='')),col=c('blue','black','gray'),
  lwd=c(3,3,3),lty=c(1,1,1),bty='n')

#### MOS vs. Aspect ####

MOS.lm.2016 <- lm(MOS1_vals[,4]~aspect_vals); s1 <- MOS.lm.2016$coefficients[2]
MOS.lm.2016a <- lm(MOS2_vals[,4]~aspect_vals); s2 <- MOS.lm.2016a$coefficients[2]
MOS.lm.2016b <- lm(MOS_vals[,4]~aspect_vals); s3 <- MOS.lm.2016b$coefficients[2]
#pdf(file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/MOS_vs_DEM.pdf')
plot(NA,xlim=c(-1,1),ylim=c(90,170),xlab='Transformed Aspect',ylab='MOS (DOY)')
abline(MOS.lm.2016,col='blue',lwd=3)
abline(MOS.lm.2016a,col='black',lwd=3)
abline(MOS.lm.2016b,col='gray',lwd=3)
legend('topleft',
  c(paste('2016 (No Corr.) (',round(s1,3),')',sep=''),
    paste('2016 (L8 Only) (',round(s2,3),')',sep=''),
    paste('2016 (HLS) (',round(s3,3),')',sep='')),col=c('blue','black','gray'),
  lwd=c(3,3,3),lty=c(1,1,1),bty='n')

#### MOA vs. Elevation ####

MOA.lm.2016 <- lm(MOA1_vals[,4]~DEM_vals); s1 <- MOA.lm.2016$coefficients[2]
MOA.lm.2016a <- lm(MOA2_vals[,4]~DEM_vals); s2 <- MOA.lm.2016a$coefficients[2]
MOA.lm.2016b <- lm(MOA_vals[,4]~DEM_vals); s3 <- MOA.lm.2016b$coefficients[2]
#pdf(file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/MOA_vs_DEM.pdf')
plot(NA,xlim=c(600,1600),ylim=c(260,330),xlab='Elevation (m)',ylab='MOA (DOY)')
abline(MOA.lm.2016,col='blue',lwd=3)
abline(MOA.lm.2016a,col='black',lwd=3)
abline(MOA.lm.2016b,col='gray',lwd=3)
legend('topleft',
  c(paste('2016 (No Corr.) (',round(s1,3),')',sep=''),
    paste('2016 (L8 Only) (',round(s2,3),')',sep=''),
    paste('2016 (HLS) (',round(s3,3),')',sep='')),col=c('blue','black','gray'),
  lwd=c(3,3,3),lty=c(1,1,1),bty='n')

#### MOA vs. Aspect ####

MOA.lm.2016 <- lm(MOA1_vals[,4]~aspect_vals); s1 <- MOA.lm.2016$coefficients[2]
MOA.lm.2016a <- lm(MOA2_vals[,4]~aspect_vals); s2 <- MOA.lm.2016a$coefficients[2]
MOA.lm.2016b <- lm(MOA_vals[,4]~aspect_vals); s3 <- MOA.lm.2016b$coefficients[2]
#pdf(file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/MOA_vs_DEM.pdf')
plot(NA,xlim=c(-1,1),ylim=c(260,330),xlab='Transformed Aspect',ylab='MOA (DOY)')
abline(MOA.lm.2016,col='blue',lwd=3)
abline(MOA.lm.2016a,col='black',lwd=3)
abline(MOA.lm.2016b,col='gray',lwd=3)
legend('topleft',
  c(paste('2016 (No Corr.) (',round(s1,3),')',sep=''),
    paste('2016 (L8 Only) (',round(s2,3),')',sep=''),
    paste('2016 (HLS) (',round(s3,3),')',sep='')),col=c('blue','black','gray'),
  lwd=c(3,3,3),lty=c(1,1,1),bty='n')

dev.off()

#### Single Pixel Plot ####
load(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/',site$name,'_bands_topocorr_method',9,sep=''))

pix <- 8673
ts1 <- red[,pix]/10000
ts2 <- nir[,pix]/10000
ts3 <- ndvi[,pix]/10000
ts4 <- evi2[,pix]/10000
df <- data.frame(yr,doy,sat,ts1,ts2,ts3,ts4)
x11()
#pdf('/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/pix_9516_valley.pdf')
par(mfrow=c(4,1),mar=c(3,3,3,3))
plot(time,ts1,pch=16,cex=1.3,main='Red')
points(time[wS],ts1[wS],pch=21,col='black',bg='yellow',cex=1.3)
legend('topright',c('S2A','L8'),pch=c(21,16),pt.bg=c('yellow','black'),bty='n')
plot(time,ts2,pch=16,cex=1.3,main='NIR')
points(time[wS],ts2[wS],pch=21,col='black',bg='yellow',cex=1.3)
plot(time,ts3,pch=16,cex=1.3,main='NDVI')
points(time[wS],ts3[wS],pch=21,col='black',bg='yellow',cex=1.3)
plot(time,ts4,pch=16,cex=1.3,main='EVI2')
points(time[wS],ts4[wS],pch=21,col='black',bg='yellow',cex=1.3)
dev.off()

source(file='/usr3/graduate/emelaas/Code/R/landsat_sentinel/v1_3/MCD12Q2C6_functions_1.R')

load(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/',site$name,'_bands_topocorr_method',9,sep=''))
pix <- 20069
ts4a <- evi2[,pix]/10000
time <- as.Date(b4[,4])
yr <- b4[,1]
doy <- b4[,2]
sat <- b4[,3]
psl1 <- PlotSeriesLandsat3(ts4a[which(yr>=2014 & yr<=2016)],time[which(yr>=2014 & yr<=2016)],pheno_pars,'','EVI2')

load(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/',site$name,'_bands_topocorr_method',1,sep=''))
ts4b <- evi2[,pix]/10000
psl2 <- PlotSeriesLandsat3(ts4b[which(yr>=2014 & yr<=2016)],time[which(yr>=2014 & yr<=2016)],pheno_pars,'','EVI2')
lines(psl1[[2]],psl1[[1]],lwd=2,col=brewer.pal(5, "Blues")[4],lty=3)
abline(v=psl1[[3]][[1]][[2]],lty=3,col='darkgreen')
abline(v=psl2[[3]][[1]][[2]],lty=1,col='darkgreen')
abline(v=psl1[[3]][[1]][[5]],lty=3,col='darkred')
abline(v=psl2[[3]][[1]][[5]],lty=1,col='darkred')
legend('topleft',c('Uncorrected','Corrected'),bty='n',cex=0.8,
  lty=c(1,3),col=c(brewer.pal(5, "Blues")[4],brewer.pal(5, "Blues")[4]))

#### In situ observations ####
load('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/random_sites/coweeta_insitu')
load('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/random_sites/phenocam_dates')
#load('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/harvard_pcam_dates')

#COW

#COW2 = 6848
#COW4 = 26714
#PhenoCam = 9098

pix <- 6848
ts4 <- evi2[,pix]/10000
ground.obs <- data.frame(as.Date(substr(strptime(paste(c(2013,2014,2015),MOS_leafon[,1]), format="%Y %j"),1,10)),
  as.Date(substr(strptime(paste(c(2013,2014,2015),MOA_color[,1]), format="%Y %j"),1,10)),
  as.Date(substr(strptime(paste(c(2013,2014,2015),MOA_leafoff[,1]), format="%Y %j"),1,10)))
PlotSeriesLandsat(ts4[which(yr<=2015)],time[which(yr<=2015)],pheno_pars,'coweeta-WS2','EVI2',ground.obs)

pix <- 26714
ts4 <- evi2[,pix]/10000
ground.obs <- data.frame(as.Date(substr(strptime(paste(c(2013,2014,2015),MOS_leafon[,2]), format="%Y %j"),1,10)),
  as.Date(substr(strptime(paste(c(2013,2014,2015),MOA_color[,2]), format="%Y %j"),1,10)),
  as.Date(substr(strptime(paste(c(2013,2014,2015),MOA_leafoff[,2]), format="%Y %j"),1,10)))
PlotSeriesLandsat(ts4[which(yr<=2015)],time[which(yr<=2015)],pheno_pars,'coweeta-WS4','EVI2',ground.obs)

pix <- 9098
ts4 <- evi2[,pix]/10000
ground.obs <- data.frame(as.Date(substr(strptime(paste(c(2013,2014,2015,2016,2017),sos50[,3]), format="%Y %j"),1,10)),
  as.Date(substr(strptime(paste(c(2013,2014,2015,2016,2017),eos50[,3]), format="%Y %j"),1,10)),
  as.Date(substr(strptime(paste(c(2013,2014,2015,2016,2017),pkr90[,3]), format="%Y %j"),1,10)))
PlotSeriesLandsat(ts4[which(yr<=2017)],time[which(yr<=2017)],pheno_pars,'coweeta-PhenoCam','EVI2',ground.obs)

#x11(h=5,w=9)
pdf(h=5,w=9,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/COW/insitu_vs_HLS.pdf')
par(mfrow=c(1,2))
plot(MOS_leafon[,1],MOS_vals[6848,1:3],xlim=c(100,150),ylim=c(100,150),pch=16,col='orange',
  ylab='HLS MOS (DOY)',xlab='In Situ MOS (DOY)',cex=1.5)
points(MOS_leafon[,2],MOS_vals[26714,1:3],pch=16,col='red',cex=1.5)
points(sos50[,3],MOS_vals[9098,1:5],pch=16,col='blue',cex=1.5)
abline(0,1,lty=2)
legend('topleft',legend=c('WS2','WS4','PhenoCam'),pch=c(16,16,16),col=c('orange','red','blue'),bty='n')

y1 <- c(MOS_vals[6848,1:3],MOS_vals[26714,1:3],MOS_vals[9098,1:5]) 
x1 <- c(MOS_leafon[,1],MOS_leafon[,2],sos50[,3])
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

plot(MOA_color[,1],MOA_vals[6848,1:3],xlim=c(260,300),ylim=c(260,300),pch=16,col='orange',
  ylab='HLS MOA (DOY)',xlab='In Situ MOA (DOY)',cex=1.5)
points(MOA_color[,2],MOA_vals[26714,1:3],pch=16,col='red',cex=1.5)
points(eos50[,3],MOA_vals[9098,1:5],pch=16,col='blue',cex=1.5)
abline(0,1,lty=2)
y1 <- c(MOA_vals[6848,1:3],MOA_vals[26714,1:3],MOA_vals[9098,1:5]) 
x1 <- c(MOA_color[,1],MOA_color[,2],eos50[,3])
lm2_1 <- lmodel2(y1~x1)
r2_1 <- lm2_1$rsquare
cortest <- cor.test(x1,y1,use='pairwise.complete.obs',method='pearson')
my.p_1 <- cortest$p.value
rmse_1 <- sqrt(mean((y1-x1)^2,na.rm=TRUE))
rp = vector('expression',2)
rp[1] = substitute(expression(paste(italic(R)^2 == MYVALUE1,' (',italic(p) == MYVALUE2,')')),
  list(MYVALUE1 = format(r2_1,dig=2),MYVALUE2 = format(my.p_1,dig=2)))[2]
rp[2] = substitute(expression(RMSE == MYVALUE),list(MYVALUE = format(rmse_1,dig=2)))[2]
legend(260,300, legend = rp, bty = 'n',cex=0.9)

dev.off()

