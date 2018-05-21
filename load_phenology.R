require(raster)
require(rgdal)
require(gdalUtils)
require(rgeos)
require(jsonlite)
require(spatialEco)
library(foreach)
library(iterators)
library(doParallel)
library(landsat)
require(data.table)
require(lubridate)

args = commandArgs(trailingOnly=T)
m <- as.numeric(args[1])
i = as.numeric(args[2])

i <- 2 #site number
m <- 5 #topo method

#List of study site locations
setwd('/projectnb/modislc/projects/landsat_sentinel/sites/')
d <- dir()
site_info <- rbind(c('COW','coweeta',35.0597,-83.4280,17,'17SKU'),
  c('HBB','hubbard',43.9438,-71.701,18,'18TYP'))

#Generate table with headers
site_info <- data.frame(site_info[,1],site_info[,2],
  as.numeric(site_info[,3]),as.numeric(site_info[,4]),as.numeric(site_info[,5]),site_info[,6])
colnames(site_info) <- c('code','name','lat','lon','UTM','granule')
site <- site_info[i,]

# Load sample image for processing
setwd(paste('/projectnb/modislc/projects/landsat_sentinel/sites/',
  site$code,'/EVI2',sep=''))
in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx("evi2*"),full.names=T,include.dirs=T,recursive=TRUE)
img <- raster(in_dirs_tile[1])

setwd(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/',site$code,sep=""))
all_files <- dir()
model <- unlist(lapply(all_files, 
  function(x) na.omit(as.numeric(unlist(strsplit(unlist(x), "[^0-9]+"))))[3]))
w <- which(model==m)
chunk <- unlist(lapply(all_files, 
  function(x) na.omit(as.numeric(unlist(strsplit(unlist(x), "[^0-9]+"))))[2]))[w]
tmp_files <- all_files[w]

phen_all <- matrix(NA,ncell(img),3)
spr2013_all <- matrix(NA,ncell(img),14)
spr2014_all <- matrix(NA,ncell(img),14)
spr2015_all <- matrix(NA,ncell(img),14)
spr2016_all <- matrix(NA,ncell(img),14)
spr2017_all <- matrix(NA,ncell(img),14)

for (j in 1:length(tmp_files)){
  print(j)
  load(tmp_files[j])
  
  phen <- matrix(unlist(pheno_mat_all[[1]]),ncol=3,byrow=TRUE)
  spr2013 <- matrix(unlist(pheno_mat_all[[2]]),ncol=14,byrow=TRUE)
  spr2014 <- matrix(unlist(pheno_mat_all[[3]]),ncol=14,byrow=TRUE)
  spr2015 <- matrix(unlist(pheno_mat_all[[4]]),ncol=14,byrow=TRUE)
  spr2016 <- matrix(unlist(pheno_mat_all[[5]]),ncol=14,byrow=TRUE)
  spr2017 <- matrix(unlist(pheno_mat_all[[6]]),ncol=14,byrow=TRUE)
  
  phen_all[(30*3660*(chunk[j]-1)+1):((30*3660)*chunk[j]),] <- phen
  spr2013_all[(30*3660*(chunk[j]-1)+1):((30*3660)*chunk[j]),] <- spr2013
  spr2014_all[(30*3660*(chunk[j]-1)+1):((30*3660)*chunk[j]),] <- spr2014
  spr2015_all[(30*3660*(chunk[j]-1)+1):((30*3660)*chunk[j]),] <- spr2015
  spr2016_all[(30*3660*(chunk[j]-1)+1):((30*3660)*chunk[j]),] <- spr2016
  spr2017_all[(30*3660*(chunk[j]-1)+1):((30*3660)*chunk[j]),] <- spr2017
}

names <- c('nobs','ncyc_gup','ncyc_gdown')
for (k in 1:3){
  print(k)
  s <- setValues(img,phen_all[,k])
  name <- paste(site_info$code[i],'_',m,'_',names[k],sep="")
  writeRaster(s,filename=paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/',
    name,sep=""),format='GTiff',overwrite=TRUE)
}

names <- c('SOS','MOS','EOS','bline', 'hline', 'gup_rsq', 'gup_avg', 'gup_longest',
  'SOA','MOA','EOA','gdown_rsq','gdown_avg','gdown_longest')
for (k in 1:14){
  print(k)
  
  if (k %in% c(1,2,3,9,10,11)){
    s1 <- setValues(img,spr2013_all[,k]-floor(unclass(as.POSIXct('2013-01-01'))/86400)[1])
    s2 <- setValues(img,spr2014_all[,k]-floor(unclass(as.POSIXct('2014-01-01'))/86400)[1])
    s3 <- setValues(img,spr2015_all[,k]-floor(unclass(as.POSIXct('2015-01-01'))/86400)[1])
    s4 <- setValues(img,spr2016_all[,k]-floor(unclass(as.POSIXct('2016-01-01'))/86400)[1])
    s5 <- setValues(img,spr2017_all[,k]-floor(unclass(as.POSIXct('2017-01-01'))/86400)[1])
  } else if (k %in% c(4,5)) {
    s1 <- setValues(img,round(spr2013_all[,k]*10000))
    s2 <- setValues(img,round(spr2014_all[,k]*10000))
    s3 <- setValues(img,round(spr2015_all[,k]*10000))
    s4 <- setValues(img,round(spr2016_all[,k]*10000))
    s5 <- setValues(img,round(spr2017_all[,k]*10000))
  } else if (k %in% c(6,12)) {
    s1 <- setValues(img,round(spr2013_all[,k]*100))
    s2 <- setValues(img,round(spr2014_all[,k]*100))
    s3 <- setValues(img,round(spr2015_all[,k]*100))
    s4 <- setValues(img,round(spr2016_all[,k]*100))
    s5 <- setValues(img,round(spr2017_all[,k]*100))
  } else {
    s1 <- setValues(img,spr2013_all[,k])
    s2 <- setValues(img,spr2014_all[,k])
    s3 <- setValues(img,spr2015_all[,k])
    s4 <- setValues(img,spr2016_all[,k])
    s5 <- setValues(img,spr2017_all[,k])
  }
  
  name <- paste(site_info$code[i],'_',m,'_',names[k],sep="")
  writeRaster(s1,filename=paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/',
    name,'_2013',sep=""),format='GTiff',overwrite=TRUE)
  writeRaster(s2,filename=paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/',
    name,'_2014',sep=""),format='GTiff',overwrite=TRUE)
  writeRaster(s3,filename=paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/',
    name,'_2015',sep=""),format='GTiff',overwrite=TRUE)
  writeRaster(s4,filename=paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/',
    name,'_2016',sep=""),format='GTiff',overwrite=TRUE)
  writeRaster(s5,filename=paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/topocorr/evi2/',
    name,'_2017',sep=""),format='GTiff',overwrite=TRUE)
}

rm(list = ls())
