library(rgdal)
library(raster)
library(foreach)
library(iterators)
library(doParallel)
library(rgeos)
require(lubridate)

all_meth <- c("No Correction","C Correction", "SCS + C", "Statistical-Empirical", "Rotation")

wshed <- readOGR('/projectnb/modislc/users/emelaas/scratch32/US-Hub/SHP/wsheds/','wsheds')
wshed_proj <- spTransform(wshed,CRS('+proj=utm +zone=18 +ellps=WGS84 +units=m +no_defs'))
wshed_buff <- gBuffer(wshed_proj,byid=FALSE,id=NULL,width=1000)

t <- read.csv('/projectnb/modislc/users/emelaas/scratch32/US-Hub/Hubbard_Brook_Data/HB_phen_coordinates.csv',
  header=TRUE)
coordinates(t) <- ~Lon+Lat
proj4string(t) <- CRS("+proj=longlat +datum=WGS84")
t_proj <- spTransform(t,CRS('+proj=utm +zone=18 +ellps=WGS84 +units=m +no_defs'))
writeOGR(t_proj,'/projectnb/modislc/projects/landsat_sentinel/shapefiles/',
  'hubbard_insitu',driver="ESRI Shapefile",overwrite=TRUE)

t_proj <- t_proj[1:9,]

for (i in 1:5){
  print(i)
  
  setwd('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/')
  in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste("HBB*",i,"*MOS*tif",sep="")),
    full.names=T,include.dirs=T,recursive=TRUE)
  MOS <- stack(in_dirs_tile)
  MOS_plots <- extract(MOS,t_proj)
  assign(paste('MOS',i,sep=""),MOS_plots)
  
  in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste("HBB*",i,"*MOA*tif",sep="")),
    full.names=T,include.dirs=T,recursive=TRUE)
  MOA <- stack(in_dirs_tile)
  MOA_plots <- extract(MOA,t_proj)
  assign(paste('MOA',i,sep=""),MOA_plots)
}

my.colors = colorRampPalette(c("red","lemonchiffon","blue"))
my.colors2 = colorRampPalette(c("brown","white","darkgreen"))
my.colors3 = colorRampPalette(c("white","gray","black"))
my.colors3b = colorRampPalette(c("black","gray","white"))
my.colors4 = colorRampPalette(c("darkgreen","yellow","orange","white"))
my.colors5 = colorRampPalette(c("white","red"))

in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste("HBB**MOS*2016.tif",sep="")),
  full.names=T,include.dirs=T,recursive=TRUE)
MOS <- stack(in_dirs_tile)
MOS_crop <- crop(MOS,wshed_buff)

in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx(paste("HBB**MOA*2016.tif",sep="")),
  full.names=T,include.dirs=T,recursive=TRUE)
MOA <- stack(in_dirs_tile)
MOA_crop <- crop(MOA,wshed_buff)

lc <- raster('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/HBB_lc.tif')
lc_crop <- crop(lc,MOS_crop[[1]])

elev <- raster('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/HBB_elev.tif')
elev_crop <- crop(elev,MOS_crop[[1]])

deg2rad <- function(deg) {(deg * pi) / (180)}
aspect <- raster('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/HBB_aspect.tif')
aspect_crop <- crop(aspect,MOS_crop[[1]])
aspect_trans <- cos(deg2rad(aspect_crop))

slope <- raster('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/HBB_slope.tif')
slope_crop <- crop(slope,MOS_crop[[1]])

nobs <- raster('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/HBB_nobs.tif')
nobs_crop <- crop(nobs,MOS_crop[[1]])

site_elev <- extract(elev_crop,t_proj)
site_aspect <- extract(aspect_trans,t_proj)


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
#pdf(h=6.5,w=12,'/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/topo_storyboard/HBB_rot_diff.pdf')
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



#x11(h=6.1,w=11)
pdf(h=6.0,w=11,'/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/topo_storyboard/HBB_maps.pdf')
par(mar=c(2,2,3.5,4),mfrow=c(2,3))
plot(elev_crop,axes=FALSE,main='Elevation (m)',legend.width=2,legend.shrink=1,
  cex.main=1.75,col=my.colors4(150),horizontal=TRUE)
plot(t_proj,add=TRUE,pch=21,col='black',bg='green',cex=1.5)

plot(aspect_trans,axes=FALSE,main='Transformed Aspect',legend.width=2,legend.shrink=1,
  cex.main=1.75,col=my.colors3(150),horizontal=TRUE)
plot(t_proj,add=TRUE,pch=21,col='black',bg='green',cex=1.5)

plot(slope_crop,axes=FALSE,main='Slope (degrees)',legend.width=2,legend.shrink=1,
  cex.main=1.75,col=my.colors3b(150),horizontal=TRUE)
plot(t_proj,add=TRUE,pch=21,col='black',bg='green',cex=1.5)

plot(nobs_crop,axes=FALSE,main='No. Observations',legend.width=2,legend.shrink=1,
  cex.main=1.75,col=my.colors5(150),horizontal=TRUE)
plot(t_proj,add=TRUE,pch=21,col='black',bg='green',cex=1.5)

plot(MOS_crop[[1]],axes=FALSE,zlim=c(120,160),colNA='black',main='Uncorrected MOS (DOY)',
  col=my.colors(150),legend.width=2,legend.shrink=1,cex.main=1.75,horizontal=TRUE)
plot(t_proj,add=TRUE,pch=21,col='black',bg='yellow',cex=1.5)

plot(MOA_crop[[1]],axes=FALSE,zlim=c(260,320),colNA='black',main='Uncorrected MOA (DOY)',
  col=my.colors(150),legend.width=2,legend.shrink=1,cex.main=1.75,horizontal=TRUE)
plot(t_proj,add=TRUE,pch=21,col='black',bg='yellow',cex=1.5)
dev.off()

#x11(h=4.7,w=12)
pdf(h=4.7,w=12,'/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/topo_storyboard/HBB_anom.pdf')
par(mar=c(2,2,4,2),mfrow=c(2,4))
for (i in 2:5){
  if (i==2) {
    plot(MOS_crop[[i]]-MOS_crop[[1]],axes=FALSE,zlim=c(-10,10),colNA='black',
      main=paste(all_meth[i],'(SPR)'),
      col=my.colors2(150),legend.width=2,legend.shrink=1,cex.main=1.75)
    plot(t_proj,add=TRUE,pch=21,col='black',bg='yellow',cex=1.5)
  } else {
    plot(MOS_crop[[i]]-MOS_crop[[1]],axes=FALSE,zlim=c(-10,10),colNA='black',
      main=paste(all_meth[i],'(SPR)'),
      col=my.colors2(150),legend=FALSE,cex.main=1.75)
    plot(t_proj,add=TRUE,pch=21,col='black',bg='yellow',cex=1.5)
  }
}
for (i in 2:5){
  if (i==2) {
    plot(MOA_crop[[i]]-MOA_crop[[1]],axes=FALSE,zlim=c(-10,10),colNA='black',
      main=paste(all_meth[i],'(AUT)'),
      col=my.colors2(150),legend=FALSE,cex.main=1.75)
    plot(t_proj,add=TRUE,pch=21,col='black',bg='yellow',cex=1.5)
  } else {
    plot(MOA_crop[[i]]-MOA_crop[[1]],axes=FALSE,zlim=c(-10,10),colNA='black',
      main=paste(all_meth[i],'(AUT)'),
      col=my.colors2(150),legend=FALSE,cex.main=1.75)
    plot(t_proj,add=TRUE,pch=21,col='black',bg='yellow',cex=1.5)
  }
}
dev.off()



#Read in Hubbard Brook ground obs.
#SPP: 1 = American Beech (FAGR); 2 = Sugar Maple (ACSA); 3 = Yellow Birch (BEAL)
setwd('/projectnb/modislc/users/emelaas/scratch32/US-Hub')
r <- read.table('HB_Phen_1989-2017.txt',header=TRUE)
r[r<0]<-NA

r_spr <- r[which(r$SEASON==1),]
r_aut <- r[which(r$SEASON==2),]

# Loop through each species, each year, and each plot and 
# estimate spring/autumn phenology dates
for (spp in 1:3){
  
  pSPR <- matrix(NA,9,5)
  pAUT <- matrix(NA,9,5)

  #Spring
  r_spr.s <- r_spr[which(r_spr$SPECIES==spp),]
  for (y in 2013:2017){
    r_spr.sy <- r_spr.s[which(r_spr.s$YEAR==y),]
    
    for (plot in 5:13){
      rint <- approx(r_spr.sy$DAY,r_spr.sy[,plot],xout=seq(1,max(r_spr.sy$DAY)))
      pSPR[(plot-4),(y-2012)] <- which.min(abs(rint$y-3))
    }
  }
  assign(paste('pSPR',spp,sep=""),pSPR)
  
  #Autumn
  r_aut.s <- r_aut[which(r_aut$SPECIES==spp),]
  for (y in 2013:2017){
    r_aut.sy <- r_aut.s[which(r_aut.s$YEAR==y),]
    
    for (plot in 5:13){
      rint <- approx(r_aut.sy$DAY,r_aut.sy[,plot],xout=seq(1,max(r_aut.sy$DAY)))
      pAUT[(plot-4),(y-2012)] <- which.min(abs(rint$y-2))
    }
  }
  assign(paste('pAUT',spp,sep=""),pAUT)
}

pSPR <- (pSPR1+pSPR2+pSPR3)/3
pAUT <- (pAUT1+pAUT2+pAUT3)/3

rmse_s1 <- sqrt(mean((as.numeric(pSPR[,3:4])-as.numeric(MOS1[1:9,3:4]))^2,na.rm=TRUE))
rmse_s2 <- sqrt(mean((as.numeric(pSPR[,3:4])-as.numeric(MOS2[1:9,3:4]))^2,na.rm=TRUE))
rmse_s3 <- sqrt(mean((as.numeric(pSPR[,3:4])-as.numeric(MOS3[1:9,3:4]))^2,na.rm=TRUE))
rmse_s4 <- sqrt(mean((as.numeric(pSPR[,3:4])-as.numeric(MOS4[1:9,3:4]))^2,na.rm=TRUE))
rmse_s5 <- sqrt(mean((as.numeric(pSPR[,3:4])-as.numeric(MOS5[1:9,3:4]))^2,na.rm=TRUE))

rmse_a1 <- sqrt(mean((as.numeric(pAUT[,3:4])-as.numeric(MOA1[1:9,3:4]))^2,na.rm=TRUE))
rmse_a2 <- sqrt(mean((as.numeric(pAUT[,3:4])-as.numeric(MOA2[1:9,3:4]))^2,na.rm=TRUE))
rmse_a3 <- sqrt(mean((as.numeric(pAUT[,3:4])-as.numeric(MOA3[1:9,3:4]))^2,na.rm=TRUE))
rmse_a4 <- sqrt(mean((as.numeric(pAUT[,3:4])-as.numeric(MOA4[1:9,3:4]))^2,na.rm=TRUE))
rmse_a5 <- sqrt(mean((as.numeric(pAUT[,3:4])-as.numeric(MOA5[1:9,3:4]))^2,na.rm=TRUE))

smat <- cbind(c(rmse_s1,rmse_s2,rmse_s3,rmse_s4,rmse_s5))
amat <- cbind(c(rmse_a1,rmse_a2,rmse_a3,rmse_a4,rmse_a5))

colors <- c('black','red','orange','yellow','white')

#x11(h=6,w=8)
pdf(h=6,w=8,'/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/topo_storyboard/HBB_barplots.pdf')
par(mfrow=c(1,2),mar=c(4,4,4,1))
barplot(smat,beside=TRUE,col=colors,ylim=c(0,10),main='Spring',
  ylab='RMSE (days)')
barplot(amat,beside=TRUE,col=colors,ylim=c(0,10),main='Autumn')
legend('topright',legend=all_meth,
  cex=0.8,fill=colors,ncol=1)
dev.off()

meanSPR <- rowMeans(cbind(pSPR1[,4],pSPR2[,4],pSPR3[,4]))
meanAUT <- rowMeans(cbind(pAUT1[,4],pAUT2[,4],pAUT3[,4]))

w <- which(getValues(lc_crop)==41)
w2 <- which(getValues(lc_crop)==41 & getValues(elev_crop)>=500 & getValues(elev_crop)<=600 & getValues(slope_crop)>=20)

#x11(h=7.2,w=7.2)
pdf(h=7.2,w=7.2,'/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/topo_storyboard/HBB_elev_aspect.pdf')
par(mfrow=c(2,2),mar=c(4,4,2,2))

MOS1.lm.2016 <- lm(getValues(MOS_crop[[1]])[w]~getValues(elev_crop)[w]); s1 <- MOS1.lm.2016$coefficients[2]
MOS2.lm.2016 <- lm(getValues(MOS_crop[[2]])[w]~getValues(elev_crop)[w]); s2 <- MOS2.lm.2016$coefficients[2]
MOS3.lm.2016 <- lm(getValues(MOS_crop[[3]])[w]~getValues(elev_crop)[w]); s3 <- MOS3.lm.2016$coefficients[2]
MOS4.lm.2016 <- lm(getValues(MOS_crop[[4]])[w]~getValues(elev_crop)[w]); s4 <- MOS4.lm.2016$coefficients[2]
MOS5.lm.2016 <- lm(getValues(MOS_crop[[5]])[w]~getValues(elev_crop)[w]); s5 <- MOS5.lm.2016$coefficients[2]
plot(NA,xlim=c(100,1600),ylim=c(130,170),xlab='Elevation (m)',ylab='MOS (DOY)')
abline(MOS1.lm.2016,col='blue',lwd=3)
abline(MOS2.lm.2016,col='green',lwd=3)
abline(MOS3.lm.2016,col='black',lwd=3)
abline(MOS4.lm.2016,col='gray',lwd=3)
abline(MOS5.lm.2016,col='orange',lwd=3)
points(site_elev[1:9],meanSPR,pch=16)
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
plot(NA,xlim=c(-1,1),ylim=c(130,150),xlab='Transformed Aspect',ylab='MOS (DOY)')
abline(MOS1.lm.2016,col='blue',lwd=3)
abline(MOS2.lm.2016,col='green',lwd=3)
abline(MOS3.lm.2016,col='black',lwd=3)
abline(MOS4.lm.2016,col='gray',lwd=3)
abline(MOS5.lm.2016,col='orange',lwd=3)
#points(site_aspect[1:9],meanSPR,pch=16)
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
plot(NA,xlim=c(100,1600),ylim=c(270,320),xlab='Elevation (m)',ylab='MOA (DOY)')
abline(MOA1.lm.2016,col='blue',lwd=3)
abline(MOA2.lm.2016,col='green',lwd=3)
abline(MOA3.lm.2016,col='black',lwd=3)
abline(MOA4.lm.2016,col='gray',lwd=3)
abline(MOA5.lm.2016,col='orange',lwd=3)
points(site_elev[1:9],meanAUT,pch=16)
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
#points(site_aspect[1:9],meanAUT,pch=16)
legend('topright',
  c(paste('No Corr (',round(s1,3),')',sep=''),
    paste('C Corr (',round(s2,3),')',sep=''),
    paste('SCS + C (',round(s3,3),')',sep=''),
    paste('S-E (',round(s4,3),')',sep=''),
    paste('Rot (',round(s5,3),')',sep='')),
  col=c('blue','green','black','gray','orange'),
  lwd=c(3,3,3,3,3),lty=c(1,1,1,1,1),bty='n')
dev.off()