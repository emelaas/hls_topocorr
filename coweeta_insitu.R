library(rgdal)
library(raster)
library(foreach)
library(iterators)
library(doParallel)
library(rgeos)
require(lubridate)

source('~/Code/R/functions/image_nan_better.R')

all_meth <- c("No Correction","C Correction", "SCS + C", "Statistical-Empirical", "Rotation")

#List of study site locations
setwd('/projectnb/modislc/projects/landsat_sentinel/sites/')
d <- dir()
site_info <- rbind(c('COW2','coweeta_sm2',35.063241,-83.432681,17,'17SKU'),
  c('COW4','coweeta_sm4',35.031436,-83.459895,17,'17SKU'),
  c('COW','coweeta_pcam',35.0597,-83.4280,17,'17SKU'))

#Generate table with headers
site_info <- data.frame(site_info[,1],site_info[,2],
  as.numeric(site_info[,3]),as.numeric(site_info[,4]),as.numeric(site_info[,5]),site_info[,6])
colnames(site_info) <- c('code','name','lat','lon','UTM','granule')

#Generate single point data frame and convert to SpatialPoint
coordinates(site_info) <- ~lon+lat
proj4string(site_info) <- CRS("+proj=longlat +datum=WGS84")
site_info_proj <- spTransform(site_info,CRS(paste('+proj=utm +zone=',site_info$UTM,
  " +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",sep="")))
writeOGR(site_info_proj,'/projectnb/modislc/projects/landsat_sentinel/shapefiles/',
  'coweeta_insitu',driver="ESRI Shapefile",overwrite=TRUE)
writeOGR(site_info_proj[3,],'/projectnb/modislc/projects/landsat_sentinel/shapefiles/',
  'coweeta_pcam',driver="ESRI Shapefile",overwrite=TRUE)

site_info_proj <- site_info_proj[1:2,]

wshed <- readOGR('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/SHP/coweeta_subwatersheds/',
  'coweeta_subwatersheds')
wshed_proj <- spTransform(wshed,CRS('+proj=utm +zone=17 +ellps=WGS84 +units=m +no_defs'))
wshed_buff <- gBuffer(wshed_proj,byid=FALSE,id=NULL,width=500)

lc <- raster('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/COW_lc.tif')
lc_crop <- crop(lc,wshed_buff)

elev <- raster('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/COW_elev.tif')
elev_crop <- crop(elev,wshed_buff)

deg2rad <- function(deg) {(deg * pi) / (180)}
aspect <- raster('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/COW_aspect.tif')
aspect_crop <- crop(aspect,wshed_buff)
aspect_trans <- cos(deg2rad(aspect_crop))

slope <- raster('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/COW_slope.tif')
slope_crop <- crop(slope,wshed_buff)

nobs <- raster('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/COW_nobs.tif')
nobs_crop <- crop(nobs,wshed_buff)

site_elev <- extract(elev_crop,site_info_proj)
site_aspect <- extract(aspect_trans,site_info_proj)


for (i in 1:5){
  print(i)
  
  setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
  in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste("COW*",i,"*MOS*tif",sep="")),
    full.names=T,include.dirs=T,recursive=TRUE)
  MOS <- stack(in_dirs_tile)
  MOS_plots <- extract(MOS,site_info_proj)
  assign(paste('MOS',i,sep=""),MOS_plots)
  
  in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste("COW*",i,"*MOA*tif",sep="")),
    full.names=T,include.dirs=T,recursive=TRUE)
  MOA <- stack(in_dirs_tile)
  MOA_plots <- extract(MOA,site_info_proj)
  assign(paste('MOA',i,sep=""),MOA_plots)
}

my.colors = colorRampPalette(c("red","lemonchiffon","blue"))
my.colors2 = colorRampPalette(c("brown","white","darkgreen"))
my.colors3 = colorRampPalette(c("white","gray","black"))
my.colors3b = colorRampPalette(c("black","gray","white"))
my.colors4 = colorRampPalette(c("darkgreen","yellow","orange","white"))
my.colors5 = colorRampPalette(c("white","red"))
my.colors6 = colorRampPalette(c("brown","yellow","darkgreen"))

in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste("COW**MOS*2016.tif",sep="")),
  full.names=T,include.dirs=T,recursive=TRUE)
MOS <- stack(in_dirs_tile)
MOS_crop <- crop(MOS,wshed_buff)

in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste("COW**MOA*2016.tif",sep="")),
  full.names=T,include.dirs=T,recursive=TRUE)
MOA <- stack(in_dirs_tile)
MOA_crop <- crop(MOA,wshed_buff)


MOS_diff <- MOS_crop[[5]]-MOS_crop[[1]]
MOS_diff_vals <- getValues(MOS_diff)
MOA_diff <- MOA_crop[[5]]-MOA_crop[[1]]
MOA_diff_vals <- getValues(MOA_diff)
ystep <- seq(-1,1,0.08)
xstep <- seq(0,50,2)
mat1 <- matrix(NA,25,25)
mat2 <- matrix(NA,25,25)
for (i in 1:25){
  for (j in 1:25){
    w <- which(getValues(aspect_trans) < ystep[i+1] & getValues(aspect_trans) >= ystep[i] & getValues(slope_crop) < xstep[j+1] & getValues(slope_crop >= xstep[j]))
    mat1[i,j] <- median(MOS_diff_vals[w],na.rm=TRUE)
    mat2[i,j] <- median(MOA_diff_vals[w],na.rm=TRUE)
  }
}
r1 <- raster(ncol=25, nrow=25, xmn=-1, xmx=1, ymn=-1, ymx=1)
r1 <- setValues(r1,mat1)
r2 <- raster(ncol=25, nrow=25, xmn=-1, xmx=1, ymn=-1, ymx=1)
r2 <- setValues(r2,mat2); r2[r2>14] <- 14
x11(h=6.5,w=12)
#pdf(h=6.5,w=12,'/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/topo_storyboard/COW_rot_diff.pdf')
par(mfrow=c(1,2))
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 4))
plot(r1,col=my.colors2(150),colNA='black',zlim=c(-14,14),axes=TRUE,
  xaxt='n',bty='n',box=FALSE,xlab='',ylab='Aspect',legend.width=1.2)
axis(1,outer=TRUE,tick=TRUE,pos=-1,las=1,at=seq(-1,1,2/5),labels=c('','10','20','30','40',''))
text(0,-1.4,'Slope')
text(0,1.3,'Difference - Middle of Spring',cex=1.3,font = 2)
plot(r2,col=my.colors2(150),colNA='black',zlim=c(-14,14),axes=TRUE,
  xaxt='n',bty='n',box=FALSE,xlab='',ylab='Aspect',legend.width=1.2)
axis(1,outer=TRUE,tick=TRUE,pos=-1,las=1,at=seq(-1,1,2/5),labels=c('','10','20','30','40',''))
text(0,-1.4,'Slope')
text(0,1.3,'Difference - Middle of Autumn',cex=1.3,font = 2)
dev.off()


#x11(h=6.8,w=9)
pdf(h=6.8,w=9,'/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/topo_storyboard/COW_maps.pdf')
par(mar=c(2,2,3.5,4),mfrow=c(2,3))
plot(elev_crop,axes=FALSE,main='Elevation (m)',legend.width=2,legend.shrink=1,
  cex.main=1.75,col=my.colors4(150),horizontal=TRUE)
plot(site_info_proj,add=TRUE,pch=21,col='black',bg='green',cex=1.5)

plot(aspect_trans,axes=FALSE,main='Transformed Aspect',legend.width=2,legend.shrink=1,
  cex.main=1.75,col=my.colors3(150),horizontal=TRUE)
plot(site_info_proj,add=TRUE,pch=21,col='black',bg='green',cex=1.5)

plot(slope_crop,axes=FALSE,main='Slope (degrees)',legend.width=2,legend.shrink=1,
  cex.main=1.75,col=my.colors3b(150),horizontal=TRUE)
plot(site_info_proj,add=TRUE,pch=21,col='black',bg='green',cex=1.5)

plot(nobs_crop,axes=FALSE,main='No. Observations',legend.width=2,legend.shrink=1,
  cex.main=1.75,col=my.colors5(150),horizontal=TRUE)
plot(site_info_proj,add=TRUE,pch=21,col='black',bg='green',cex=1.5)

plot(MOS_crop[[1]],axes=FALSE,zlim=c(90,150),colNA='black',main='Uncorrected MOS (DOY)',
  col=my.colors(150),legend.width=2,legend.shrink=1,cex.main=1.75,horizontal=TRUE)
plot(site_info_proj,add=TRUE,pch=21,col='black',bg='yellow',cex=1.5)

plot(MOA_crop[[1]],axes=FALSE,zlim=c(260,330),colNA='black',main='Uncorrected MOA (DOY)',
  col=my.colors(150),legend.width=2,legend.shrink=1,cex.main=1.75,horizontal=TRUE)
plot(site_info_proj,add=TRUE,pch=21,col='black',bg='yellow',cex=1.5)
dev.off()

#x11(h=6,w=11)
pdf(h=5.9,w=11,'/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/topo_storyboard/COW_anom.pdf')
par(mar=c(2,2,4,2),mfrow=c(2,4))
for (i in 2:5){
  if (i==2) {
    plot(MOS_crop[[i]]-MOS_crop[[1]],axes=FALSE,zlim=c(-10,10),colNA='black',
      main=paste(all_meth[i],'(SPR)'),
      col=my.colors2(150),legend.width=2,legend.shrink=1,cex.main=1.75)
    plot(site_info_proj,add=TRUE,pch=21,col='black',bg='yellow',cex=1.5)
  } else {
    plot(MOS_crop[[i]]-MOS_crop[[1]],axes=FALSE,zlim=c(-10,10),colNA='black',
      main=paste(all_meth[i],'(SPR)'),
      col=my.colors2(150),legend=FALSE,cex.main=1.75)
    plot(site_info_proj,add=TRUE,pch=21,col='black',bg='yellow',cex=1.5)
  }
}
for (i in 2:5){
  if (i==2) {
    plot(MOA_crop[[i]]-MOA_crop[[1]],axes=FALSE,zlim=c(-10,10),colNA='black',
      main=paste(all_meth[i],'(AUT)'),
      col=my.colors2(150),legend=FALSE,cex.main=1.75)
    plot(site_info_proj,add=TRUE,pch=21,col='black',bg='yellow',cex=1.5)
  } else {
    plot(MOA_crop[[i]]-MOA_crop[[1]],axes=FALSE,zlim=c(-10,10),colNA='black',
      main=paste(all_meth[i],'(AUT)'),
      col=my.colors2(150),legend=FALSE,cex.main=1.75)
    plot(site_info_proj,add=TRUE,pch=21,col='black',bg='yellow',cex=1.5)
  }
}
dev.off()


w <- which(getValues(lc_crop)==41)
w2 <- which(getValues(lc_crop)==41 & getValues(elev_crop)>=900 & getValues(elev_crop)<=1000 & getValues(slope_crop)>=20)

#x11(h=7.2,w=7.2)
pdf(h=7.2,w=7.2,'/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/topo_storyboard/COW_elev_aspect.pdf')
par(mfrow=c(2,2),mar=c(4,4,2,2))

MOS1.lm.2016 <- lm(getValues(MOS_crop[[1]])[w]~getValues(elev_crop)[w]); s1 <- MOS1.lm.2016$coefficients[2]
MOS2.lm.2016 <- lm(getValues(MOS_crop[[2]])[w]~getValues(elev_crop)[w]); s2 <- MOS2.lm.2016$coefficients[2]
MOS3.lm.2016 <- lm(getValues(MOS_crop[[3]])[w]~getValues(elev_crop)[w]); s3 <- MOS3.lm.2016$coefficients[2]
MOS4.lm.2016 <- lm(getValues(MOS_crop[[4]])[w]~getValues(elev_crop)[w]); s4 <- MOS4.lm.2016$coefficients[2]
MOS5.lm.2016 <- lm(getValues(MOS_crop[[5]])[w]~getValues(elev_crop)[w]); s5 <- MOS5.lm.2016$coefficients[2]
plot(NA,xlim=c(600,1600),ylim=c(100,150),xlab='Elevation (m)',ylab='MOS (DOY)')
abline(MOS1.lm.2016,col='blue',lwd=3)
abline(MOS2.lm.2016,col='green',lwd=3)
abline(MOS3.lm.2016,col='black',lwd=3)
abline(MOS4.lm.2016,col='gray',lwd=3)
abline(MOS5.lm.2016,col='orange',lwd=3)
legend('topleft',
  c(paste('No Corr (',round(s1,3),')',sep=''),
    paste('C Corr (',round(s2,3),')',sep=''),
    paste('SCS + C (',round(s3,3),')',sep=''),
    paste('S-E (',round(s4,3),')',sep=''),
    paste('Rot (',round(s5,3),')',sep='')),
    col=c('blue','green','black','gray','orange'),
  lwd=c(3,3,3,3,3),lty=c(1,1,1,1,1),bty='n')

MOS1.lm.2016 <- lm(getValues(MOS_crop[[1]])[w2]~getValues(aspect_trans)[w2]); s1 <- MOS1.lm.2016$coefficients[2]
MOS2.lm.2016 <- lm(getValues(MOS_crop[[2]])[w2]~getValues(aspect_trans)[w2]); s2 <- MOS2.lm.2016$coefficients[2]
MOS3.lm.2016 <- lm(getValues(MOS_crop[[3]])[w2]~getValues(aspect_trans)[w2]); s3 <- MOS3.lm.2016$coefficients[2]
MOS4.lm.2016 <- lm(getValues(MOS_crop[[4]])[w2]~getValues(aspect_trans)[w2]); s4 <- MOS4.lm.2016$coefficients[2]
MOS5.lm.2016 <- lm(getValues(MOS_crop[[5]])[w2]~getValues(aspect_trans)[w2]); s5 <- MOS5.lm.2016$coefficients[2]
plot(NA,xlim=c(-1,1),ylim=c(100,140),xlab='Transformed Aspect',ylab='MOS (DOY)')
abline(MOS1.lm.2016,col='blue',lwd=3)
abline(MOS2.lm.2016,col='green',lwd=3)
abline(MOS3.lm.2016,col='black',lwd=3)
abline(MOS4.lm.2016,col='gray',lwd=3)
abline(MOS5.lm.2016,col='orange',lwd=3)
legend('topleft',
  c(paste('No Corr (',round(s1,3),')',sep=''),
    paste('C Corr (',round(s2,3),')',sep=''),
    paste('SCS + C (',round(s3,3),')',sep=''),
    paste('S-E (',round(s4,3),')',sep=''),
    paste('Rot (',round(s5,3),')',sep='')),
  col=c('blue','green','black','gray','orange'),
  lwd=c(3,3,3,3,3),lty=c(1,1,1,1,1),bty='n')

MOA1.lm.2016 <- lm(getValues(MOA_crop[[1]])[w]~getValues(elev_crop)[w]); s1 <- MOA1.lm.2016$coefficients[2]
MOA2.lm.2016 <- lm(getValues(MOA_crop[[2]])[w]~getValues(elev_crop)[w]); s2 <- MOA2.lm.2016$coefficients[2]
MOA3.lm.2016 <- lm(getValues(MOA_crop[[3]])[w]~getValues(elev_crop)[w]); s3 <- MOA3.lm.2016$coefficients[2]
MOA4.lm.2016 <- lm(getValues(MOA_crop[[4]])[w]~getValues(elev_crop)[w]); s4 <- MOA4.lm.2016$coefficients[2]
MOA5.lm.2016 <- lm(getValues(MOA_crop[[5]])[w]~getValues(elev_crop)[w]); s5 <- MOA5.lm.2016$coefficients[2]
plot(NA,xlim=c(600,1600),ylim=c(280,310),xlab='Elevation (m)',ylab='MOA (DOY)')
abline(MOA1.lm.2016,col='blue',lwd=3)
abline(MOA2.lm.2016,col='green',lwd=3)
abline(MOA3.lm.2016,col='black',lwd=3)
abline(MOA4.lm.2016,col='gray',lwd=3)
abline(MOA5.lm.2016,col='orange',lwd=3)
legend('topleft',
  c(paste('No Corr (',round(s1,3),')',sep=''),
    paste('C Corr (',round(s2,3),')',sep=''),
    paste('SCS + C (',round(s3,3),')',sep=''),
    paste('S-E (',round(s4,3),')',sep=''),
    paste('Rot (',round(s5,3),')',sep='')),
  col=c('blue','green','black','gray','orange'),
  lwd=c(3,3,3,3,3),lty=c(1,1,1,1,1),bty='n')

MOA1.lm.2016 <- lm(getValues(MOA_crop[[1]])[w2]~getValues(aspect_trans)[w2]); s1 <- MOA1.lm.2016$coefficients[2]
MOA2.lm.2016 <- lm(getValues(MOA_crop[[2]])[w2]~getValues(aspect_trans)[w2]); s2 <- MOA2.lm.2016$coefficients[2]
MOA3.lm.2016 <- lm(getValues(MOA_crop[[3]])[w2]~getValues(aspect_trans)[w2]); s3 <- MOA3.lm.2016$coefficients[2]
MOA4.lm.2016 <- lm(getValues(MOA_crop[[4]])[w2]~getValues(aspect_trans)[w2]); s4 <- MOA4.lm.2016$coefficients[2]
MOA5.lm.2016 <- lm(getValues(MOA_crop[[5]])[w2]~getValues(aspect_trans)[w2]); s5 <- MOA5.lm.2016$coefficients[2]
plot(NA,xlim=c(-1,1),ylim=c(280,310),xlab='Transformed Aspect',ylab='MOA (DOY)')
abline(MOA1.lm.2016,col='blue',lwd=3)
abline(MOA2.lm.2016,col='green',lwd=3)
abline(MOA3.lm.2016,col='black',lwd=3)
abline(MOA4.lm.2016,col='gray',lwd=3)
abline(MOA5.lm.2016,col='orange',lwd=3)
legend('topright',cex=1,
  c(paste('No Corr (',round(s1,3),')',sep=''),
    paste('C Corr (',round(s2,3),')',sep=''),
    paste('SCS + C (',round(s3,3),')',sep=''),
    paste('S-E (',round(s4,3),')',sep=''),
    paste('Rot (',round(s5,3),')',sep='')),
  col=c('blue','green','black','gray','orange'),
  lwd=c(3,3,3,3,3),lty=c(1,1,1,1,1),bty='n')
dev.off()






load('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/random_sites/coweeta_insitu')
load('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/phenocam_dates')



#### SPRING ####

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/csv/COW/')
dat <- read.table('Spring_1142_1_1142.CSV',skip=5,sep=',')
colnames(dat) <- c('site','branch_no','branch_ID','height','species','date','bud_stage','comments')
dat$date <- as.Date(strptime(dat$date,"%m/%d/%Y"))

yr <- as.numeric(substring(dat$date,1,4))
dat <- dat[which(yr>=2013),]

yr <- as.numeric(substr(dat$date,1,4))
mo <- as.numeric(substr(dat$date,6,7))
day <- as.numeric(substr(dat$date,9,10))
doy <- yday(as.Date(paste(yr,'-',mo,'-',day,sep="")))

dat$yr <- yr
dat$doy <- doy

MOS_leafon <- matrix(NA,3,2)
for (i in c(2,4)){
  dat_site <- dat[which(dat$site==i),]
  
  for (yr in c(2013:2015)){
    dat_site_yr <- dat_site[which(dat_site$yr==yr),]
    x <- dat_site_yr$doy
    y <- dat_site_yr$bud_stage
    m <- nls(y ~ 1+4/(1+exp(-c*(x-d))), start = list(c=0.1,d=120))
    m.summary <- summary(m)
    m.coef <- m.summary$coefficients[,1]
    m.pred <- 1+4/(1+exp(-m.coef[1]*(seq(1,150)-m.coef[2])))
    MOS_leafon[yr-2012,i/2] <- which.min(abs(m.pred-4))
  }
}

#### AUTUMN ####

setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/csv/COW/')
dat <- read.table('Fall_1142_1_1142.CSV',skip=5,sep=',')
colnames(dat) <- c('site','branch_no','branch_ID','height','species','date','color','leafoff','comments')
dat$date <- as.Date(strptime(dat$date,"%m/%d/%Y"))

yr <- as.numeric(substring(dat$date,1,4))
dat <- dat[which(yr>=2013),]

yr <- as.numeric(substr(dat$date,1,4))
mo <- as.numeric(substr(dat$date,6,7))
day <- as.numeric(substr(dat$date,9,10))
doy <- yday(as.Date(paste(yr,'-',mo,'-',day,sep="")))

dat$yr <- yr
dat$doy <- doy

MOA_color <- matrix(NA,3,2)
MOA_leafoff <- matrix(NA,3,2)
for (i in c(2,4)){
  dat_site <- dat[which(dat$site==i),]
  
  for (yr in c(2013:2015)){
    dat_site_yr <- dat_site[which(dat_site$yr==yr),]
    x <- dat_site_yr$doy
    y <- dat_site_yr$color
    m <- nls(y ~ 100/(1+exp(-c*(x-d))), start = list(c=0.1,d=260))
    m.summary <- summary(m)
    m.coef <- m.summary$coefficients[,1]
    m.pred <- 100/(1+exp(-m.coef[1]*(seq(1,350)-m.coef[2])))
    MOA_color[yr-2012,i/2] <- which.min(abs(m.pred-50))
    
    x <- dat_site_yr$doy
    y <- dat_site_yr$leafoff
    m <- nls(y ~ 100/(1+exp(-c*(x-d))), start = list(c=0.1,d=260))
    m.summary <- summary(m)
    m.coef <- m.summary$coefficients[,1]
    m.pred <- 100/(1+exp(-m.coef[1]*(seq(1,350)-m.coef[2])))
    MOA_leafoff[yr-2012,i/2] <- which.min(abs(m.pred-50))
  }
}

save(MOA_leafoff,MOA_color,MOS_leafon,file='/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/coweeta_insitu')


rmse_s11 <- sqrt(mean((as.numeric(t(MOS_leafon))-as.numeric(MOS1[1:2,1:3]))^2,na.rm=TRUE))
rmse_s12 <- sqrt(mean((as.numeric(t(MOS_leafon))-as.numeric(MOS2[1:2,1:3]))^2,na.rm=TRUE))
rmse_s13 <- sqrt(mean((as.numeric(t(MOS_leafon))-as.numeric(MOS3[1:2,1:3]))^2,na.rm=TRUE))
rmse_s14 <- sqrt(mean((as.numeric(t(MOS_leafon))-as.numeric(MOS4[1:2,1:3]))^2,na.rm=TRUE))
rmse_s15 <- sqrt(mean((as.numeric(t(MOS_leafon))-as.numeric(MOS5[1:2,1:3]))^2,na.rm=TRUE))

rmse_c11 <- sqrt(mean((as.numeric(t(MOA_color))-as.numeric(MOA1[1:2,1:3]))^2,na.rm=TRUE))
rmse_c12 <- sqrt(mean((as.numeric(t(MOA_color))-as.numeric(MOA2[1:2,1:3]))^2,na.rm=TRUE))
rmse_c13 <- sqrt(mean((as.numeric(t(MOA_color))-as.numeric(MOA3[1:2,1:3]))^2,na.rm=TRUE))
rmse_c14 <- sqrt(mean((as.numeric(t(MOA_color))-as.numeric(MOA4[1:2,1:3]))^2,na.rm=TRUE))
rmse_c15 <- sqrt(mean((as.numeric(t(MOA_color))-as.numeric(MOA5[1:2,1:3]))^2,na.rm=TRUE))

rmse_a11 <- sqrt(mean((as.numeric(t(MOA_leafoff))-as.numeric(MOA1[1:2,1:3]))^2,na.rm=TRUE))
rmse_a12 <- sqrt(mean((as.numeric(t(MOA_leafoff))-as.numeric(MOA2[1:2,1:3]))^2,na.rm=TRUE))
rmse_a13 <- sqrt(mean((as.numeric(t(MOA_leafoff))-as.numeric(MOA3[1:2,1:3]))^2,na.rm=TRUE))
rmse_a14 <- sqrt(mean((as.numeric(t(MOA_leafoff))-as.numeric(MOA4[1:2,1:3]))^2,na.rm=TRUE))
rmse_a15 <- sqrt(mean((as.numeric(t(MOA_leafoff))-as.numeric(MOA5[1:2,1:3]))^2,na.rm=TRUE))


mat <- cbind(c(rmse_s11,rmse_s12,rmse_s13,rmse_s14,rmse_s15),
  c(rmse_c11,rmse_c12,rmse_c13,rmse_c14,rmse_c15),
  c(rmse_a11,rmse_a12,rmse_a13,rmse_a14,rmse_a15))

colors <- c('black','red','orange','yellow','white')

pdf('/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/topo_storyboard/COW_barplots.pdf')
barplot(mat,beside=TRUE,col=colors,ylim=c(0,15),names=c('Leaf On','Leaf Color','Leaf Off'),
  ylab='RMSE (days)')
legend('topleft',legend=all_meth,
  cex=0.8,fill=colors,ncol=1)
dev.off()