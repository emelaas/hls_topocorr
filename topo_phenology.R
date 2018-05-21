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

#Register the parallel backend
registerDoParallel(8)

args = commandArgs(trailingOnly=T)
method <- as.numeric(args[1])
chunk = as.numeric(args[2])

print(c(method,chunk))

sgranule <- readOGR('/projectnb/modislc/projects/landsat_sentinel/shapefiles/sentinel2_tiles_world/sentinel2_tiles_world/','sentinel2_tiles_world')

#Choose site
for (i in 1:2){
  print(i)
        
  comb <- function(x, ...) {
    lapply(seq_along(x),
      function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
    
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
  
  source(file='/usr3/graduate/emelaas/Code/GitHub/hls/MCD12Q2C6_functions_1.R')
  pheno_pars <- list(
    LandsatFillQuant=0.05,
    LandsatXmin=0,
    LandsatSpikeThresh=2,
    LandsatMinResid=0.1,
    LandsatFillDOY=NULL,
    LandsatDoAnnual=T,
    LandsatPadHeadTail=T,
    min_peak_to_peak_distance=50,
    min_peak_quantile=0.2,
    max_seg_length=200,
    min_seg_amplitude=0.00,
    agg_amp_frac=0.15,
    gup_threshes=c(0.1,0.5,0.9),
    gdown_threshes=c(0.9,0.5,0.1),
    spline_spar=0.5
  )
    
  #Find all images for each site
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/sites/',
    site$code,'/EVI2/',method,sep=''))
  in_dirs_tile <- list.files(path=getwd(),pattern=glob2rx("evi2*tif"),full.names=T,include.dirs=T,recursive=TRUE)
  
  r <- stack(in_dirs_tile)
  c <- crop(r,extent(r, 30*(chunk-1)+1, 30*(chunk), 1, 3660))
  vi_vals <- getValues(c)
  vi_vals[vi_vals<=0] <- NA
  
  LC <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site_info$code[i],'_lc.tif',sep=''))
  LC_crop <- crop(LC,extent(r, 30*(chunk-1)+1, 30*(chunk), 1, 3660))
  LC_vals <- getValues(LC_crop)
  
  #Generate index of ancillary info for all images
  doy <- as.numeric(substring(in_dirs_tile,74,76))
  yr <- as.numeric(substring(in_dirs_tile,70,73))
  sat <- substring(in_dirs_tile,68,68)
  
  index <- data.frame(sat,yr,doy,t(vi_vals))
  colnames(index) <- c('sat','yr','doy','vi')
  index <- index[order(yr,doy),]
  
  vi_vals <- index[,-c(1,2,3)]
  info <- index[,1:3]
  
  # Convert to data table to speed up VI time series extraction (below)
  dt.evi2 <- data.table(vi_vals,keep.rownames=FALSE)
  
  pheno_mat_all <- foreach(i = 1:ncol(vi_vals), .combine='comb', .multicombine=TRUE,
    .init=list(list(), list(), list(), list(), list(), list())) %dopar% {
      if (i%%10000==0) print(i)
      
      all <- list(as.numeric(matrix(NA,14,1)),as.numeric(matrix(NA,14,1)),
        as.numeric(matrix(NA,14,1)),as.numeric(matrix(NA,14,1)),as.numeric(matrix(NA,14,1)))
      
      evi2_ts <- as.matrix(dt.evi2[,..i]/10000)
      index2 <- cbind(info,evi2_ts)
      time <- as.Date(substr(strptime(paste(info[,2], info[,3]), format="%Y %j"),1,10))
      daynum <- as.numeric(time-time[1])
      #index3 <- index2[-which(is.na(index2[,4])==1 | index2[,3]=='S'),]
      #time <- time[-which(is.na(index2[,4])==1 | index2[,3]=='S')]
      
      pheno_metrics <- try(DoPhenologyLandsat(as.matrix(index2[,4]),time,pheno_pars))
      
      if (length(which(is.na(index2[,4])==0))>0 & is.na(pheno_metrics)==0){
        num_obs <- length(which(is.na(index2[,4])==0)) #Calculate number of cloud-free observations
        w_gup <- which(pheno_metrics[[8]]-pheno_metrics[[7]]>0.1) #which green up shoulders have amplitude > 0.10?
        w_gdown <- which(pheno_metrics[[10]]-pheno_metrics[[9]]>0.1) #which green down shoulders?
        
        num_cyc_gup <- length(w_gup) #total number of green ups
        num_cyc_gdown <- length(w_gdown) #total number of green downs
        
        if (min(c(num_cyc_gup,num_cyc_gdown))>0){
          all_yrs <- as.numeric(substr(pheno_metrics[[1]],1,4))
          w_gup_yrs <- all_yrs[w_gup]-2012
          w_gdown_yrs <- all_yrs[w_gdown]-2012
          
          
          # spr1, spr2, spr3, seg_min_gup, seg_max_gup, gup_rsq, gup_missing, gup_longest
          tmp <- unlist(pheno_metrics)
          for (j in 1:length(w_gup)){
            a <- tmp[seq(w_gup[j],length(tmp),length(all_yrs))] #green up metrics
            all[[w_gup_yrs[j]]][c(1:8)] <- a[c(1,2,3,7,8,12,13,14)]
          }
          
          # aut1, aut2, aut3, gdown_rsq, gdown_missing, gdown_longest
          for (j in 1:length(w_gdown)){
            b <- tmp[seq(w_gdown[j],length(tmp),length(all_yrs))] #green down metrics
            all[[w_gdown_yrs[j]]][c(9:14)] <- b[c(4,5,6,15,16,17)]
          }
        } else {
          num_obs <- length(which(is.na(index2[,4])==0))
        }
        
      } else {
        num_obs <- length(which(is.na(index2[,4])==0))
        num_cyc_gup <- 0
        num_cyc_gdown <- 0
      }
      
      list(c(num_obs,num_cyc_gup,num_cyc_gdown),all[[1]],all[[2]],all[[3]],all[[4]],all[[5]])
    }
    
  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/',site$code,sep=""))
  print(paste('Saving ',site$code,'_','evi2_topocorr_phenology_',chunk,'_',method,sep=""))
  save(pheno_mat_all,
    file=paste(site$code,'_','evi2_topocorr_phenology_',chunk,'_',method,sep=""))  
}
