require(rgdal)
library(raster)

stack_loc <- 'ah_hubbard'

#Read in inventory data (Schwarz et al. 2003 Ecology)
setwd('/projectnb/modislc/users/emelaas/scratch32/US-Hub/')
Lon <- read.csv('LON.csv',header=FALSE)
Lat <- read.csv('LAT.csv',header=FALSE)
crown_status <- read.csv('CROWN_STATUS.csv',header=FALSE)
plot <- read.csv('PLOT.csv',header=FALSE)
dbh <- read.csv('DBH.csv',header=FALSE)
spp <- read.table('spp.csv',header=FALSE)

data <- data.frame(plot,Lat,Lon,spp,dbh,crown_status)
colnames(data) <- c('plot','Lat','Lon','spp','dbh','crown_status')

#Sugar Maple, Yellow Birch, American Beech
#Count number of Dominant, Codominant, and Intermediante stems of each species 
#and make a table ('spp_perc')
w_ALL <- which(data[,1] >= 1)
w_CD <- which(data[,6]=='C' | data[,6]=='D' | data[,6]=='I')
w_SM <-which(data[,4]=='SM' & data[,6]=='C' | data[,4]=='SM' & data[,6]=='D' | data[,4]=='SM' & data[,6]=='I')
w_YB <-which(data[,4]=='YB' & data[,6]=='C' | data[,4]=='YB' & data[,6]=='D' | data[,4]=='YB' & data[,6]=='I')
w_AB <-which(data[,4]=='AB' & data[,6]=='C' | data[,4]=='AB' & data[,6]=='D' | data[,4]=='AB' & data[,6]=='I')

spp_perc <- matrix(0,460,8)
spp_perc[as.numeric(rownames(table(data[w_ALL,1]))),1] <- as.matrix(table(data[w_ALL,1]))
spp_perc[as.numeric(rownames(table(data[w_CD,1]))),2] <- as.matrix(table(data[w_CD,1]))
spp_perc[as.numeric(rownames(table(data[w_SM,1]))),3] <- as.matrix(table(data[w_SM,1]))
spp_perc[as.numeric(rownames(table(data[w_YB,1]))),4] <- as.matrix(table(data[w_YB,1]))
spp_perc[as.numeric(rownames(table(data[w_AB,1]))),5] <- as.matrix(table(data[w_AB,1]))
spp_perc[,6] <- spp_perc[,3]+spp_perc[,4]+spp_perc[,5]

lat <- aggregate(data$Lat,by=list(data$plot),FUN=median)
lon <- aggregate(data$Lon,by=list(data$plot),FUN=median)
spp_perc[lat[,1],7] <- lat[,2]
spp_perc[lon[,1],8] <- lon[,2]
spp_perc <- spp_perc[-(which(spp_perc[,1]==0)),]

inventory <- data.frame(spp_perc[,1],spp_perc[,7],spp_perc[,8])
colnames(inventory) <- c('plot','lat','lon')

#Load in LC and DEM Map
setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,sep=""))
if (file.exists('EOSD')==1){
  lc <- raster(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/EOSD/eosd_clip.bip',sep=""))
} else if (file.exists('NLCD')==1) {
  lc <- raster(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/NLCD/nlcd_clip.bip',sep=""))
}

setwd('/projectnb/modislc/users/emelaas/scratch32/US-Hub/USGS-DEM/')
DEM <- raster('NE_DEM_UTM_reproj.tif')
DEM_crop <- crop(DEM,lc)

#Generate Pixel Number map using LC map
N <- seq(1,ncell(lc),1)
tmp <- setValues(lc,N)

#Generate table of pheno station characteristics for extracting dates 
sta <- c('1B','6T','4B','4T','5B','5T','7B','7T')
lat <- c(43.952022,43.955547,43.950917,43.958947,43.949128,43.957678,43.928928,43.918358)
lon <- c(-71.725558,-71.742178,-71.728739,-71.731661,-71.731831,-71.736869,-71.765647,-71.769833)
elev <- c(472,767,528,608,505,727,605,823)
plots <- data.frame(sta,lat,lon,elev)

#Create Spatial Data Frame for phenology plots
pheno_plots <- plots
coordinates(pheno_plots) <- ~lon+lat
proj4string(pheno_plots) <- CRS("+proj=longlat +datum=WGS84")
pheno_plots_proj <- spTransform(pheno_plots,CRS("+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs"))

pheno_extr <- data.frame(unlist(extract(tmp,pheno_plots_proj)),seq(1,8))
colnames(pheno_extr) <- c('pixID','plot')

#Convert Lat/Lon of inventory pixels to UTM coordinates
coordinates(inventory) <- ~lon+lat
proj4string(inventory) <- CRS("+proj=longlat +datum=WGS84")
inv_proj <- spTransform(inventory,CRS("+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs"))

inv_extr <- data.frame(unlist(extract(tmp,inv_proj)))
inv_DEM_extr <- data.frame(unlist(extract(DEM_crop,inv_proj)))
inv_extr <- data.frame(inv_extr,inv_DEM_extr,spp_perc)
colnames(inv_extr) <- c('pixID','DEM','ALL','CD','SM','YB','AB','SUM','Lat','Lon')

#Load in AH Hubbard Landsat phenology matrix
setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/phenology',sep=""))
load('pheno_nodist_mat')
pheno_mat <- data.frame(pheno_mat)
names(pheno_mat)[1] <- 'pixID'

#Merge phenology_matrix with pixels co-located with phenology plots
pheno_pix <- merge(pheno_mat,pheno_extr,by='pixID')
pheno_pix[pheno_pix==0] <- NA
pheno_pix <- pheno_pix[order(pheno_pix$plot),]

setwd('/projectnb/modislc/users/emelaas/scratch32/US-Hub')
save(pheno_pix,file='pheno_pix')

pheno_spr <- as.matrix(pheno_pix[,17:40])
pheno_aut <- as.matrix(pheno_pix[,49:72])

#Merge phenology matrix with pixels co-located with inventory plots
inv_pix <- merge(pheno_mat,inv_extr,by='pixID')
inv_pix[inv_pix==0] <- NA

pix_spr <- inv_pix[,17:40]
pix_aut <- inv_pix[,49:72]


#Read in Hubbard Brook ground obs.
#SPP: 1 = American Beech (FAGR); 2 = Sugar Maple (ACSA); 3 = Yellow Birch (BEAL)
setwd('/projectnb/modislc/users/emelaas/scratch32/US-Hub')
r <- read.table('HB_Phen_1989-2017.txt',header=TRUE)
r[r<0]<-NA

r_spr <- r[which(r$SEASON==1),]
r_aut <- r[which(r$SEASON==2),]

#Loop through each species, each year, and each plot and 
# 1. estimate spring/autumn phenology dates
# 2. calculate RMSE, MBE, and COR between ground obs. and satellite pixel at each plot
for (spp in 1:3){
  
  pSPR <- matrix(NA,8,29)
  pAUT <- matrix(NA,8,29)
#   RMSEspr <- matrix(NA,428,8)
#   RMSEaut <- matrix(NA,428,8)
#   MBEspr <- matrix(NA,428,8)
#   MBEaut <- matrix(NA,428,8)
#   CORspr <- matrix(NA,428,8)
#   CORaut <- matrix(NA,428,8)
  
  #Spring
  r_spr.s <- r_spr[which(r_spr$SPECIES==spp),]
  for (y in 1989:2017){
      r_spr.sy <- r_spr.s[which(r_spr.s$YEAR==y),]
      
    for (plot in 5:12){
      rint <- approx(r_spr.sy$DAY,r_spr.sy[,plot],xout=seq(1,max(r_spr.sy$DAY)))
      pSPR[(plot-4),(y-1988)] <- which.min(abs(rint$y-3))
    }
  }
  assign(paste('pSPR',spp,sep=""),pSPR)
  
#   #Calculate Spring statistics
#   for (inv in 1:428){
#     for (plot in 1:8){
#       obs <- as.numeric(pSPR[plot,])
#       sat <- as.numeric(pix_spr[inv,])
#       RMSEspr[inv,plot] <- sqrt(mean((obs-sat)^2,na.rm=TRUE))
#       MBEspr[inv,plot] <- mean(obs-sat,na.rm=TRUE)
#       CORspr[inv,plot] <- cor(obs,sat,use='pairwise.complete.obs')
#     }
#   }
#   assign(paste('RMSEspr',spp,sep=""),RMSEspr)
#   assign(paste('MBEspr',spp,sep=""),MBEspr)
#   assign(paste('CORspr',spp,sep=""),CORspr)
  
  #Autumn
  r_aut.s <- r_aut[which(r_aut$SPECIES==spp),]
  for (y in 1989:2017){
    r_aut.sy <- r_aut.s[which(r_aut.s$YEAR==y),]
    
    for (plot in 5:12){
      rint <- approx(r_aut.sy$DAY,r_aut.sy[,plot],xout=seq(1,max(r_aut.sy$DAY)))
      pAUT[(plot-4),(y-1988)] <- which.min(abs(rint$y-2))
    }
  }
  assign(paste('pAUT',spp,sep=""),pAUT)
  
#   #Calculate autumn statistics
#   for (inv in 1:428){
#     for (plot in 1:8){
#       obs <- as.numeric(pAUT[plot,])
#       sat <- as.numeric(pix_aut[inv,])
#       RMSEaut[inv,plot] <- sqrt(mean((obs-sat)^2,na.rm=TRUE))
#       MBEaut[inv,plot] <- mean(obs-sat,na.rm=TRUE)
#       CORaut[inv,plot] <- cor(obs,sat,use='pairwise.complete.obs') 
#     }
#   }
#   assign(paste('RMSEaut',spp,sep=""),RMSEaut)
#   assign(paste('MBEaut',spp,sep=""),MBEaut)
#   assign(paste('CORaut',spp,sep=""),CORaut)
}

#Find pixels dominated by each species (more than 75% of all C/D/I)
w1 <- which(inv_pix$AB/inv_pix$CD>0.75)
w2 <- which(inv_pix$SM/inv_pix$CD>0.75)
w3 <- which(inv_pix$YB/inv_pix$CD>0.75)

#Calculate distance between each phenostation and each inventory cell
dist <- matrix(NA,428,8)
for (i in 1:428){
  for (j in 1:8){
    dist[i,j] <- sqrt((plots$lat[j]-inv_pix$Lat[i])^2+(plots$lon[j]-inv_pix$Lon[i])^2)
  }
}
w <- apply(dist,2,which.min)
sta_pix <- inv_pix[w,]

#Compare species-weighted ground phenology (based on sta_pix above) with aggregate
#Landsat phenology for corresponding MODIS pixel

setwd(paste('/projectnb/modislc/projects/te_phenology/landsat_stacks/',stack_loc,'/phenology',sep=""))
load('panel_nodist_stats')
panel_stats <- data.frame(panel_stats)
names(panel_stats)[1] <- 'modisID'
names(sta_pix)[4] <- 'modisID'

inv_panel <- merge(panel_stats,sta_pix,by='modisID')
panel_spr <- as.matrix(inv_panel[,13:36])
panel_aut <- as.matrix(inv_panel[,45:68])

sta_pix[is.na(sta_pix)==1] <- 0
sta_spr <- matrix(NA,8,24)
sta_aut <- matrix(NA,8,24)
for (i in 1:8){
  sta_spr[i,] <- (sta_pix$AB/sta_pix$SUM)*pSPR1[i,] + (sta_pix$SM/sta_pix$SUM)*pSPR2[i,] + 
    (sta_pix$YB/sta_pix$SUM)*pSPR3[i,]
  sta_aut[i,] <- (sta_pix$AB/sta_pix$SUM)*pAUT1[i,] + (sta_pix$SM/sta_pix$SUM)*pAUT2[i,] + 
    (sta_pix$YB/sta_pix$SUM)*pAUT3[i,]
}

panel_spr_all <- panel_spr
dim(panel_spr_all) <- c(24*8,1)
panel_aut_all <- panel_aut
dim(panel_aut_all) <- c(24*8,1)
sta_spr_all <- sta_spr
dim(sta_spr_all) <- c(24*8,1)
sta_aut_all <- sta_aut
dim(sta_aut_all) <- c(24*8,1)

stats_spr_fagr <- matrix(NA,8,3)
stats_spr_acsa <- matrix(NA,8,3)
stats_spr_beal <- matrix(NA,8,3)
stats_aut_fagr <- matrix(NA,8,3)
stats_aut_acsa <- matrix(NA,8,3)
stats_aut_beal <- matrix(NA,8,3)

for (i in 1:8) {
  stats_spr_fagr[i,1] <- sqrt(mean((pheno_spr[i,]-pSPR1[i,])^2,na.rm=TRUE))
  stats_spr_fagr[i,2] <- mean(pheno_spr[i,]-pSPR1[i,],na.rm=TRUE) #Landsat minus ground
  stats_spr_fagr[i,3] <- cor(pheno_spr[i,],pSPR1[i,],use='pairwise.complete.obs')^2
  stats_spr_acsa[i,1] <- sqrt(mean((pheno_spr[i,]-pSPR2[i,])^2,na.rm=TRUE))
  stats_spr_acsa[i,2] <- mean(pheno_spr[i,]-pSPR2[i,],na.rm=TRUE)
  stats_spr_acsa[i,3] <- cor(pheno_spr[i,],pSPR2[i,],use='pairwise.complete.obs')^2
  stats_spr_beal[i,1] <- sqrt(mean((pheno_spr[i,]-pSPR3[i,])^2,na.rm=TRUE))
  stats_spr_beal[i,2] <- mean(pheno_spr[i,]-pSPR3[i,],na.rm=TRUE)
  stats_spr_beal[i,3] <- cor(pheno_spr[i,],pSPR3[i,],use='pairwise.complete.obs')^2
  
  stats_aut_fagr[i,1] <- sqrt(mean((pheno_aut[i,]-pAUT1[i,])^2,na.rm=TRUE))
  stats_aut_fagr[i,2] <- mean(pheno_aut[i,]-pAUT1[i,],na.rm=TRUE)
  stats_aut_fagr[i,3] <- cor(pheno_aut[i,],pAUT1[i,],use='pairwise.complete.obs')^2
  stats_aut_acsa[i,1] <- sqrt(mean((pheno_aut[i,]-pAUT2[i,])^2,na.rm=TRUE))
  stats_aut_acsa[i,2] <- mean(pheno_aut[i,]-pAUT2[i,],na.rm=TRUE)
  stats_aut_acsa[i,3] <- cor(pheno_aut[i,],pAUT2[i,],use='pairwise.complete.obs')^2
  stats_aut_beal[i,1] <- sqrt(mean((pheno_aut[i,]-pAUT3[i,])^2,na.rm=TRUE))
  stats_aut_beal[i,2] <- mean(pheno_aut[i,]-pAUT3[i,],na.rm=TRUE)
  stats_aut_beal[i,3] <- cor(pheno_aut[i,],pAUT3[i,],use='pairwise.complete.obs')^2
}

#x11(h=5.5,w=8)
pdf(h=6,w=8,'/projectnb/modislc/projects/te_phenology/figures/Storyboard/ah_hubbard_sat_v_ground_boxplots.pdf')
par(mfrow=c(2,3),mar=c(3,4.5,2.25,0.25))
boxplot(stats_spr_acsa[,1],stats_spr_beal[,1],stats_spr_fagr[,1],names=c('ACSA','BEAL','FAGR'),col=c('red','white','blue'),
  ylab='RMSE (days)',ylim=c(0,18))
mtext('(a)',side=3,line=0.5,cex=0.85,adj=-0.1,font=2)
boxplot(stats_spr_acsa[,2],stats_spr_beal[,2],stats_spr_fagr[,2],names=c('ACSA','BEAL','FAGR'),col=c('red','white','blue'),
  ylab='MBE (days)',main='Spring',ylim=c(-18,18))
abline(h=0,lty=2)
mtext('(b)',side=3,line=0.5,cex=0.85,adj=-0.1,font=2)
boxplot(stats_spr_acsa[,3],stats_spr_beal[,3],stats_spr_fagr[,3],names=c('ACSA','BEAL','FAGR'),col=c('red','white','blue'),
  ylab=expression(R^2),ylim=c(0,1))
mtext('(c)',side=3,line=0.5,cex=0.85,adj=-0.1,font=2)
boxplot(stats_aut_acsa[,1],stats_aut_beal[,1],stats_aut_fagr[,1],names=c('ACSA','BEAL','FAGR'),col=c('red','white','blue'),
  ylab='RMSE (days)',ylim=c(0,18))
mtext('(d)',side=3,line=0.5,cex=0.85,adj=-0.1,font=2)
boxplot(stats_aut_acsa[,2],stats_aut_beal[,2],stats_aut_fagr[,2],names=c('ACSA','BEAL','FAGR'),col=c('red','white','blue'),
  ylab='MBE (days)',main='Autumn',ylim=c(-18,18))
abline(h=0,lty=2)
mtext('(e)',side=3,line=0.5,cex=0.85,adj=-0.1,font=2)
boxplot(stats_aut_acsa[,3],stats_aut_beal[,3],stats_aut_fagr[,3],names=c('ACSA','BEAL','FAGR'),col=c('red','white','blue'),
  ylab=expression(R^2),ylim=c(0,1))
mtext('(f)',side=3,line=0.5,cex=0.85,adj=-0.1,font=2)

dev.off()


x11(h=8.5,w=8.5)
#pdf(h=8.5,w=8.5,'/projectnb/modislc/projects/te_phenology/figures/Storyboard/ah_hubbard_sat_v_ground.pdf')
par(mfrow=c(4,4),mar=c(4,4,1.25,0.25))
for (i in 1:8){
  if (i==1){
    plot(pheno_spr[i,],pSPR1[i,],pch=16,cex=1.5,col='red',xlim=c(120,165),ylim=c(120,165),
      xlab='Landsat SOS (DOY)',ylab='Ground SOS (DOY)',main='Spring')
  } else {
    plot(pheno_spr[i,],pSPR1[i,],pch=16,cex=1.5,col='red',xlim=c(120,165),ylim=c(120,165),
      xlab='',ylab='')
  }
  points(pheno_spr[i,],pSPR2[i,],pch=16,cex=1.5,col='blue',xlim=c(120,165),ylim=c(120,165))
  points(pheno_spr[i,],pSPR3[i,],pch=16,cex=1.5,col='green',xlim=c(120,165),ylim=c(120,165))
  abline(0,1,lty=2)
  mtext(paste(' ',plots$sta[i]),side=3,line=-1.5,cex=1,adj=0.01,font=2)
  
  if (i==4){
    legend('bottomright',pch=c(16,16,16),pt.cex=c(1.5,1.5,1.5),bty='n',cex=1,
      col=c('red','green','blue'),legend=c('FAGR','ACSA','BEAL'))
  }
}
for (i in 1:8){
  if (i==1){
    plot(pheno_aut[i,],pAUT1[i,],pch=16,cex=1.5,col='red',xlim=c(245,300),ylim=c(245,300),
      xlab='Landsat EOS (DOY)',ylab='Ground EOS (DOY)',main='Autumn')
  } else {
    plot(pheno_aut[i,],pAUT1[i,],pch=16,cex=1.5,col='red',xlim=c(245,300),ylim=c(245,300),
      xlab='',ylab='')
  }
  points(pheno_aut[i,],pAUT2[i,],pch=16,cex=1.5,col='blue',xlim=c(245,300),ylim=c(245,300))
  points(pheno_aut[i,],pAUT3[i,],pch=16,cex=1.5,col='green',xlim=c(245,300),ylim=c(245,300))
  abline(0,1,lty=2)
  mtext(paste(' ',plots$sta[i]),side=3,line=-1.5,cex=1,adj=0.01,font=2)
}

dev.off()

##PLOT MODIS PIXEL-AGGREGATE PHENOLOGY VS. GROUND OBS WEIGHTED BY SPP COMPOSITION OF 
##NEAREST INVENTORY PLOT (SCHWARZ DATASET)

#x11(h=4,w=8)
pdf(h=4,w=8,'/projectnb/modislc/projects/te_phenology/figures/Storyboard/ah_hubbard_sat_v_ground.pdf')
par(mfrow=c(1,2))
plot(panel_spr[1,]-mean(panel_spr[1,],na.rm=TRUE),sta_spr[1,]-mean(sta_spr[1,],na.rm=TRUE),
  pch=16,xlim=c(-20,20),ylim=c(-20,20),xlab='Landsat SOS Anomaly (Days)',ylab='Ground SOS Anomaly (Days)')
abline(0,1,lty=2)
for (i in 2:8){
  points(panel_spr[i,]-mean(panel_spr[i,],na.rm=TRUE),sta_spr[i,]-mean(sta_spr[i,],na.rm=TRUE),
    pch=16,xlim=c(125,165),ylim=c(125,165))
}
mtext('(a)',side=3,line=0.5,cex=1,adj=-0.1,font=2)

lm.spr <- lm(sta_spr_all~panel_spr_all)
lm.spr_coef <- round(coef(lm.spr),3)
r2 <- summary(lm.spr)$adj.r.squared
abline(lm.spr)
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
  list(MYVALUE = format(r2,dig=2)))[2]
# rp[2] = substitute(expression(y == MYVALUE2*x+MYVALUE3), 
#   list(MYVALUE2 = format(lm.spr_coef[2], digits = 4),
#     MYVALUE3 = format(lm.spr_coef[1], digits = 4)))[2]
legend('bottomright', legend = rp, bty = 'n')

plot(panel_aut[1,]-mean(panel_aut[1,],na.rm=TRUE),sta_aut[1,]-mean(sta_aut[1,],na.rm=TRUE),
  pch=16,xlim=c(-20,20),ylim=c(-20,20),xlab='Landsat EOS Anomaly (Days)',ylab='Ground EOS Anomaly (Days)')
abline(0,1,lty=2)
for (i in 2:8){
  points(panel_aut[i,]-mean(panel_aut[i,],na.rm=TRUE),sta_aut[i,]-mean(sta_aut[i,],na.rm=TRUE),
    pch=16,xlim=c(125,165),ylim=c(125,165))
}
mtext('(b)',side=3,line=0.5,cex=1,adj=-0.1,font=2)

lm.aut <- lm(sta_aut_all~panel_aut_all)
lm.aut_coef <- round(coef(lm.aut),3)
r2 <- summary(lm.aut)$adj.r.squared
abline(lm.aut)
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
  list(MYVALUE = format(r2,dig=2)))[2]
# rp[2] = substitute(expression(y == MYVALUE2*x+MYVALUE3), 
#   list(MYVALUE2 = format(lm.aut_coef[2], digits = 4),
#     MYVALUE3 = format(lm.aut_coef[1], digits = 4)))[2]
legend('bottomright', legend = rp, bty = 'n')

dev.off()
