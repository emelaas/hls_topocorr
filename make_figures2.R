library(lmodel2)
library(raster)
library(rgdal)

deg2rad <- function(deg) {(deg * pi) / (180)}

i <- 1
m <- 2
yr <- 2016

setwd('/projectnb/modislc/projects/landsat_sentinel/sites/')
d <- dir()
site_info <- rbind(c('COW','coweeta',35.0597,-83.4280,17,'17SKU'),
  c('HBB','hubbard',43.9438,-71.701,18,'18TYP'))

#Generate table with headers
site_info <- data.frame(site_info[,1],site_info[,2],
  as.numeric(site_info[,3]),as.numeric(site_info[,4]),as.numeric(site_info[,5]),site_info[,6])
colnames(site_info) <- c('code','name','lat','lon','UTM','granule')

site <- site_info[i,]

wshed <- raster('/projectnb/modislc/projects/landsat_sentinel/orthophoto/COW/conif_perc.tif')

LC <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site_info$code[i],'_lc.tif',sep=''))
slope <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site_info$code[i],'_slope.tif',sep=''))
aspect <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site_info$code[i],'_aspect.tif',sep=''))
aspect_trans <- cos(deg2rad(aspect))
elev <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site_info$code[i],'_DEM.tif',sep=''))

decid_pix <- which(getValues(LC)==41 & getValues(slope) >= 20)

SOS <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/',
  site_info$code[i],'_',m,'_SOS_',yr,'.tif',sep=''))
MOS <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/',
  site_info$code[i],'_',m,'_MOS_',yr,'.tif',sep=''))
EOS <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/',
  site_info$code[i],'_',m,'_EOS_',yr,'.tif',sep=''))
SOA <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/',
  site_info$code[i],'_',m,'_SOA_',yr,'.tif',sep=''))
MOA <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/',
  site_info$code[i],'_',m,'_MOA_',yr,'.tif',sep=''))
EOA <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/',
  site_info$code[i],'_',m,'_EOA_',yr,'.tif',sep=''))

round_elev <- 10*round(getValues(elev)[decid_pix]/10)
round_aspect <- round(getValues(aspect_trans)[decid_pix],1)

boxplot(getValues(MOA)[decid_pix]~round_elev,outline=FALSE)
boxplot(getValues(MOA)[decid_pix]~round_aspect,outline=FALSE)