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
registerDoParallel(16)

#Choose site
for (i in 1:2){
  print(i)
  
  args = commandArgs(trailingOnly=T)
  #m = as.numeric(args[1])
  m <- 4

  all_meth <- c("none","ccorrection", "scsc", "empirical", "rotational", "illumination")

  source(file='/usr3/graduate/emelaas/Code/GitHub/hls/topocorr_v4.R')

  comb <- function(x, ...) {
    lapply(seq_along(x),
      function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }

  sgranule <- readOGR('/projectnb/modislc/projects/landsat_sentinel/shapefiles/sentinel2_tiles_world/sentinel2_tiles_world/','sentinel2_tiles_world')

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

  setwd(paste('/projectnb/modislc/projects/landsat_sentinel/sites/',site$code,'/EVI2/',sep=''))
  dir.create(paste(getwd(),'/',m,sep=''))
  
  #Find all images for each site
  setwd('/projectnb/modislc/projects/landsat_sentinel/sites/')
  data_loc <- paste(getwd(),'/',site_info$code[i],'/bands_v1_3',sep='')
  in_dirs <- list.files(path=data_loc,pattern=glob2rx("HLS**3.tif"),full.names=T,include.dirs=T,recursive=TRUE)
  granules <- unique(substr(in_dirs,78,80))
  in_dirs_tile <- list.files(path=data_loc,pattern=glob2rx(paste("HLS*",granules[1],"*3.tif",sep="")),
    full.names=T,include.dirs=T,recursive=TRUE)

  img <- raster(in_dirs_tile[1])

  #Generate DEM, LC, slope and aspect maps
  sgranule_subset1 <- sgranule[which(sgranule$Name==levels(droplevels(site_info$granule[i]))),]
  sgranule_proj <- spTransform(sgranule_subset1,CRS(as.character(crs(img))))
  sgranule_proj <- crop(sgranule_proj,img)
  writeOGR(sgranule_proj,'/projectnb/modislc/projects/landsat_sentinel/shapefiles/',
    site_info$code[i],driver="ESRI Shapefile",overwrite=TRUE)

  setwd('/projectnb/modislc/projects/landsat_sentinel/shapefiles/')

#   src_data2 <- '/projectnb/modislc/data/lc_database/regional/united_states/NLCD2006_landcover_4-20-11_se5/nlcd2006_landcover_4-20-11_se5.img'
#   LC <- gdalwarp(src_data2,dstfile=paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site_info$code[i],'_lc.tif',sep=""),
#     t_srs=projection(sgranule_proj),ts=c(3660,3660),
#     cutline=paste(site_info$code[i],'.shp',sep=""),cl=site_info$code[i],crop_to_cutline=TRUE,output_Raster=TRUE,
#     overwrite=FALSE,verbose=TRUE)
  LC <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site_info$code[i],'_lc.tif',sep=''))
  LC_vals <- getValues(LC)
  LC_types <- as.numeric(rownames(table(LC_vals)))
  
#   src_data <- '/projectnb/modislc/data/dem/usgs_ned/mosaic_slope_aea_proj.tif'
#   slope <- gdalwarp(src_data,
#     dstfile=paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site_info$code[i],'_slope.tif',sep=""),
#     t_srs=projection(sgranule_proj),ts=c(3660,3660),
#     cutline=paste(site_info$code[i],'.shp',sep=""),cl=site_info$code[i],
#     crop_to_cutline=TRUE,output_Raster=TRUE,overwrite=FALSE,verbose=TRUE)
  slope <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site_info$code[i],'_slope.tif',sep=''))
  slope_vals <- getValues(slope)
  
  src_data <- '/projectnb/modislc/data/dem/usgs_ned/mosaic_dem_aea_proj.tif'
  elev <- gdalwarp(src_data,
    dstfile=paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site_info$code[i],'_elev.tif',sep=""),
    t_srs=projection(sgranule_proj),ts=c(3660,3660),
    cutline=paste(site_info$code[i],'.shp',sep=""),cl=site_info$code[i],
    crop_to_cutline=TRUE,output_Raster=TRUE,overwrite=FALSE,verbose=TRUE)
  elev <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site_info$code[i],'_elev.tif',sep=''))
  elev_vals <- getValues(elev)
  
#   src_data <- '/projectnb/modislc/data/dem/usgs_ned/mosaic_aspect_aea_proj.tif'
#   aspect <- gdalwarp(src_data,
#     dstfile=paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site_info$code[i],'_aspect.tif',sep=""),
#     t_srs=projection(sgranule_proj),ts=c(3660,3660),
#     cutline=paste(site_info$code[i],'.shp',sep=""),cl=site_info$code[i],
#     crop_to_cutline=TRUE,output_Raster=TRUE,overwrite=FALSE,verbose=TRUE)
  aspect <- raster(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/tifs/',site_info$code[i],'_aspect.tif',sep=''))
  aspect_vals <- getValues(aspect)

  #Compile Solar Zenith Angle time series
  data_loc2 <- '/projectnb/modislc/projects/landsat_sentinel/v1_3/csv/'
  sza_file <- list.files(path=data_loc2,pattern=glob2rx(paste(site_info$code[i],'*csv',sep="")))
  tmp <- as.numeric(unlist(read.csv(paste(data_loc2,sza_file,sep=''),header=FALSE)))
  tmp1 <- t(matrix(tmp,ncol=length(tmp)/4,byrow=TRUE))
  tmp2 <- tmp1[order(tmp1[,1],tmp1[,2]),]
  sza <- tmp2[,3]
  saa <- tmp2[,4]

  #Loop through each image in the site stack
  all_info <- foreach(j = 1:length(in_dirs_tile), .combine = rbind) %dopar% {
      
      print(j)
    
      #Crop each image using arbitrarily defined window
      r <- stack(in_dirs_tile[j])

      #Save cropped raster
      tmp <- r[[1]]
      tmp <- setValues(tmp,seq(1,ncell(slope)))

      band_vals <- getValues(r)
      band_vals[band_vals<=0] <- NA

      w <- which(band_vals[,6]!=0 & band_vals[,6]!=64 & band_vals[,6]!=128 & band_vals[,6]!=192)
      band_vals[w,1:5] <- NA
      
      b1_corr <- as.numeric(matrix(NA,ncell(tmp),1))
      b2_corr <- as.numeric(matrix(NA,ncell(tmp),1))
      b3_corr <- as.numeric(matrix(NA,ncell(tmp),1))
      b4_corr <- as.numeric(matrix(NA,ncell(tmp),1))
      b5_corr <- as.numeric(matrix(NA,ncell(tmp),1))
      
      for (n in 1:length(LC_types)){
        print(n)
        wLC <- which(LC_vals==LC_types[n] & is.na(aspect_vals)==0 & is.na(band_vals[,1])==0)
        if (length(wLC) > 4){
          b1_corr_LC <- topocorr_v4(x=band_vals[wLC,1],slope=slope_vals[wLC], 
            aspect=aspect_vals[wLC], sunelev=90-sza[j], sunazimuth=saa[j],
            method=all_meth[m])
          b1_corr[wLC] <- b1_corr_LC
        }
        
        wLC <- which(LC_vals==LC_types[n] & is.na(aspect_vals)==0 & is.na(band_vals[,2])==0)
        if (length(wLC) > 4){
          b2_corr_LC <- topocorr_v4(x=band_vals[wLC,2],slope=slope_vals[wLC], 
            aspect=aspect_vals[wLC], sunelev=90-sza[j], sunazimuth=saa[j],
            method=all_meth[m])
          b2_corr[wLC] <- b2_corr_LC
        }
        
        wLC <- which(LC_vals==LC_types[n] & is.na(aspect_vals)==0 & is.na(band_vals[,3])==0)
        if (length(wLC) > 4){
          b3_corr_LC <- topocorr_v4(x=band_vals[wLC,3],slope=slope_vals[wLC], 
            aspect=aspect_vals[wLC], sunelev=90-sza[j], sunazimuth=saa[j],
            method=all_meth[m])
          b3_corr[wLC] <- b3_corr_LC
        }
        
        wLC <- which(LC_vals==LC_types[n] & is.na(aspect_vals)==0 & is.na(band_vals[,4])==0)
        if (length(wLC) > 4){
          b4_corr_LC <- topocorr_v4(x=band_vals[wLC,4],slope=slope_vals[wLC], 
            aspect=aspect_vals[wLC], sunelev=90-sza[j], sunazimuth=saa[j],
            method=all_meth[m])
          b4_corr[wLC] <- b4_corr_LC
        }
        
        wLC <- which(LC_vals==LC_types[n] & is.na(aspect_vals)==0 & is.na(band_vals[,5])==0)
        if (length(wLC) > 4){
          b5_corr_LC <- topocorr_v4(x=band_vals[wLC,5],slope=slope_vals[wLC], 
            aspect=aspect_vals[wLC], sunelev=90-sza[j], sunazimuth=saa[j],
            method=all_meth[m])
          b5_corr[wLC] <- b5_corr_LC
        }
      }
    
      ndvi <- (b4_corr/10000 - b3_corr/10000)/(b4_corr/10000 + b3_corr/10000)
      evi2 <- 2.5*(b4_corr/10000 - b3_corr/10000)/(b4_corr/10000 + 2.4*b3_corr/10000 + 1)
    
      ndvi <- round(ndvi*10000)
      evi2 <- round(evi2*10000)
          
      #Generate index of ancillary info for all images
      time <- as.numeric(substring(in_dirs_tile[[j]],82,88))
      doy <- as.numeric(substring(in_dirs_tile[[j]],86,88))
      yr <- as.numeric(substring(in_dirs_tile[[j]],82,85))
      sat <- substring(in_dirs_tile[[j]],71,71)

      evi2 <- setValues(tmp,evi2)
    
      writeRaster(evi2,
      filename=paste('/projectnb/modislc/projects/landsat_sentinel/sites/',
        site$code,'/EVI2/',m,'/evi2_',sat,'_',time,sep=""),format='GTiff',overwrite=TRUE)
    
      info <- c(doy,yr,sat)
    }

  save(all_info,file=paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/',
    site$name,'_all_info',m,sep=''))
}
