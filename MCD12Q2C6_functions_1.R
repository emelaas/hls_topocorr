library(rgdal)
library(raster)
library(foreach)
library(iterators)
library(doParallel)
library(zoo)
library(RColorBrewer)
require(data.table)
#library(GISTools)

#args = commandArgs(trailingOnly=T)
#chunk = as.numeric(args[1]) 

# chunk <- 1

#Register the parallel backend
registerDoParallel(16)

comb <- function(x, ...) {
  lapply(seq_along(x),
    function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

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
  min_seg_amplitude=0.05,
  agg_amp_frac=0.15,
  gup_threshes=c(0.1,0.5,0.9),
  gdown_threshes=c(0.9,0.5,0.1),
  spline_spar=0.5
)

# Jamie Olson's local min/max method from here: http://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
GetLocalMaxMin <- function(x, partial=TRUE, minima=FALSE){
  if(minima){
    if(partial){
      which(diff(c(FALSE,diff(x)>0,TRUE))>0)
    }else{
      which(diff(diff(x)>0)>0)+1
    }
  }else{
    if(partial){
      which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
    }else{
      which(diff(diff(x)>=0)<0)+1
    }
  }
}
#---------------------------------------------------------------------
# Returns indices in x of potential peak values using the GetLocalMaxMin function
FindPotentialPeaks_minmax <- function(x){
  # check for all NA
  if(all(is.na(x))) return(NA)
  # check for constant x
  if((min(x, na.rm=T) - max(x, na.rm=T)) == 0) return(NA)
  potential_peaks <- try(GetLocalMaxMin(x), silent=T)
  if(inherits(potential_peaks, 'try-error')){
    return(NA)
  }else{
    # return the list of potential peaks
    return(potential_peaks)
  }
}
#---------------------------------------------------------------------
# screens peaks to make sure that they are separated by at least min_peak_to_peak_distance
# if they're not, the lesser of the two peaks is eliminated.
# valid peaks must also be greater than the min_peak_quantile of x
FilterPeaks <- function(x, potential_peaks, min_peak_to_peak_distance, min_peak_quantile=0.2){
  # check that potential_peaks is not NA
  if(all(is.na(potential_peaks))) return(NA)
  # first eliminate all peaks that have an absolute magnitude below min_peak_quantile
  # potential_peaks <- potential_peaks[!(x[potential_peaks] < quantile(x[potential_peaks], min_peak_quantile))] # among the peaks
  potential_peaks <- potential_peaks[!(x[potential_peaks] < quantile(x, min_peak_quantile))] # among the EVI2 time series
  # we loop through all the potential peaks, only moving on when we have satisfied
  # the minimum distance between peaks requirement
  i <- 1
  while(i < length(potential_peaks)){
    # find the distance to next peak
    peak_to_peak_dist <- potential_peaks[i + 1] - potential_peaks[i]
    if(peak_to_peak_dist < min_peak_to_peak_distance){
      # eliminate the smaller of the two potential peaks
      if(x[potential_peaks[i + 1]] >= x[potential_peaks[i]]){
        potential_peaks <- potential_peaks[-1 * i]
      }else{
        potential_peaks <- potential_peaks[-1 * (i + 1)]
      }
    }else{
      # distance is fine, move on
      i <- i + 1
    }
    # if we've eliminated the last peak, return NA
    if(length(potential_peaks) == 0) return(NA)
  }
  return(potential_peaks)
}
#---------------------------------------------------------------------
IsSpike <- function(x, thresh=2, minResid=0){
  # this function detects negative spikes with an ad-hoc, 3-point method.
  # the deviation of the middle point from the line between the first and last
  # point is calculated. If it is negative, the absolute value of the
  # ratio of this deviation and the change between the first and last points is
  # calculated. If the ratio is greater than "thresh", and the deviation is greater
  # than minResid, then the point is considered a negative spike and T is returned
  # y <- x
  # x <- 1:3
  # dev <- y[2] - ((((y[3] - y[1]) * (x[2] - x[1])) / (x[3] - x[1])) + y[1])
  dev <- x[2] - ((x[3] - x[1]) / 2) - x[1]
  # devRatio <- dev / (y[3] - y[1])
  devRatio <- dev / (x[3] - x[1])
  ifelse((abs(dev) > minResid) & (dev < 0) & (abs(devRatio) > thresh), T, F)
}
#---------------------------------------------------------------------
CheckSpike <- function(x, thresh=2, minResid=0){
  # applies the IsSpike function to a time series
  x_og <- x # preserve original vector
  x_outs <- rep(F, length(x_og)) # create the outlier output vector
  x <- x_og[!is.na(x_og)] # subset to non missing values
  # apply the IsSpike function on a rolling 3-point window
  outs <- rollapply(x, 3, IsSpike, thresh=thresh, minResid=minResid, partial=T)
  outs[is.na(outs)] <- F # replace NAs with FALSE
  x_outs[!is.na(x_og)] <- outs # expand to size of original x, accounting for missing values
  return(x_outs)
}
#---------------------------------------------------------------------
LandsatSmoother <- function(x, dates, fill_quant=0.05, x_min=0, spike_thresh=2, min_resid=0.1, spline_spar=0.25, fill_doy=NULL, doAnnual=T, padHeadTail=T){
  # special function to screen/smooth/fill the Landsat time series
  # if there are less than 4 valid values, we can't fit a spline
  if(sum(!is.na(x)) < 4) return(NA)
  # create daily dates and VI vector
  if(doAnnual){
    # if doAnnual is TRUE, then the pred dates are for full calendar years, regardless of the min/max of dates
    pred_dates <- seq(as.Date(paste(strftime(min(dates), format="%Y"), "-1-1", sep="")), as.Date(paste(strftime(max(dates), format="%Y"), "-12-31", sep="")), by="day")
  }else{
    pred_dates <- seq(min(dates), max(dates), by="day")
  }
  x_doy <- rep(NA, length(pred_dates))
  x_doy[as.numeric(dates)-as.numeric(pred_dates[1])+1] <- x # assign original data to proper dates
  #x_doy[pred_dates %in% dates] <- x 
  
  # despike
  x_tmp <- x_doy[!is.na(x_doy)] # remove missing values for despiking
  outliers <- CheckSpike(x_tmp, thresh=spike_thresh, minResid=min_resid)
  x_tmp[outliers] <- NA
  x_doy[!is.na(x_doy)] <- x_tmp # replace the original values with outlier screened values
  
  # fill min values
  #minVI <- quantile(x_doy[x_doy > x_min], fill_quant, na.rm=T)
  #x_doy[x_doy < minVI] <- minVI
  
  # fill dormant period if requested, each period is filled with 3 values only
  num_fill_values <- 3
  if(!is.null(fill_doy)){
    doy <- as.numeric(strftime(pred_dates, format="%j"))
    for(fill_seg in fill_doy){
      fill_doys <- c(fill_seg[1], sample((fill_seg[1] + 1):(fill_seg[2] - 1), (num_fill_values - 2)), fill_seg[2])
      x_doy[doy %in% fill_doys] <- minVI
    }
  }
  
  # experiment with weighting high values the most
  # xe <- ecdf(x_doy)
  # w <- xe(x_doy)
  # w[w < 0.7 & !is.na(w)] <- 1
  # w[w >= 0.7 & !is.na(w)] <- w[w >= 0.7 & !is.na(w)] * 10
  
  # check once more for too few non-missing values
  if(sum(!is.na(x_doy)) < 4) return(NA)
  
  # smooth with a spline to get continuous daily series
  # pred_dates <- seq(min(dates), max(dates), by="day")
  spl <- smooth.spline(pred_dates[!is.na(x_doy)], x_doy[!is.na(x_doy)], spar=spline_spar)
  # weighted version
  # spl <- smooth.spline(pred_dates[!is.na(x_doy)], x_doy[!is.na(x_doy)], w=w[!is.na(x_doy)], spar=spline_spar)
  x_smooth <- predict(spl, as.numeric(pred_dates))$y
  # screen again for low values
  #x_smooth[x_smooth < minVI] <- minVI
  
  # pad the head and tail with minVI, if requested
  if(padHeadTail){
    #x_doy[1:min(which(!is.na(x_doy)))] <- minVI
    #x_doy[max(which(!is.na(x_doy))):length(x_doy)] <- minVI
    
    x_smooth[1:min(which(!is.na(x_doy)))] <- x_doy[min(which(!is.na(x_doy)))]
    x_smooth[max(which(!is.na(x_doy))):length(x_doy)] <- x_doy[max(which(!is.na(x_doy)))]
  }
  
  return(list(pred_dates, x_smooth))
}
#---------------------------------------------------------------------
# returns ascending segments in x
GetFullSegs <- function(x, peaks, max_seg_length, min_seg_amplitude, agg_amp_frac=0.15){
  # check if peaks is NA and return NA if so
  if(all(is.na(peaks))) return(NA)
  # find ascending nodes by locating the minimum/trough value before the
  # peak and within search window, checking the amplitude
  full_segs <- NULL
  tmp_seg <- NA
  for(i in 1:length(peaks)){
    #################################
    # First, find ascending segments
    # define search window
    search_window_start <- max(
      1,
      peaks[i - 1],
      peaks[i] - max_seg_length
    )
    # find troughs in search window, return one closest to peak that meets min_seg_amplitude requirement,
    # if there are none, find minimum value
    troughs <- GetLocalMaxMin(x[search_window_start:peaks[i]], minima=T)
    if(length(troughs) == 0){
      # no troughs, use min value in search window
      potential_trough <- search_window_start + which(x[search_window_start:peaks[i]] == min(x[search_window_start:peaks[i]], na.rm=T))[1] - 1
      # check amplitude, add segment to list if it covers minimum amplitude. Otherwise, do nothing and move on
      tmp_amp <- x[peaks[i]] - x[potential_trough]
      # if(x[peaks[i]] - x[potential_trough] > min_seg_amplitude){
      if(!is.na(tmp_amp)){
        if(x[peaks[i]] - x[potential_trough] > min_seg_amplitude){
          tmp_seg <- list(c(potential_trough, peaks[i]))
        }
      }
    }else{
      # so there are troughs, we want to find the one furthest from the peak that yields a segment that is both
      # greater than min_seg_amplitude, and increases the amplitude more than the agg_amp_frac.
      seg_amp <- 0
      true_trough <- NULL
      for(trough in rev(troughs)){
        potential_trough <- trough + search_window_start
        tmp_amp <- x[peaks[i]] - x[potential_trough]
        if(!is.na(tmp_amp)){
          if((tmp_amp >= min_seg_amplitude) & (tmp_amp >= ((1 + agg_amp_frac) * seg_amp))){
            true_trough <- potential_trough
            seg_amp <- tmp_amp
          }
        }
      }
      if(!is.null(true_trough)){
        tmp_seg <- c(true_trough, peaks[i])
      }
    }# end check for troughs in search window
    #################################
    # Find end of segment
    # if we found an ascending segment, then we need to find an associated trough/min point after the peak
    if(!all(is.na(tmp_seg))){
      search_window_end <- min(
        length(x),
        peaks[i + 1], # this can result in NA when it reads off the end
        peaks[i] + max_seg_length,
        na.rm=T
      )
      # find troughs in search window, return one closest to peak that meets min_seg_amplitude requirement,
      # if there are none, find minimum value
      troughs <- GetLocalMaxMin(x[peaks[i]:search_window_end], minima=T)
      if(length(troughs) == 0){
        # no troughs, use min value in search window
        potential_trough <- peaks[i] + which(x[peaks[i]:search_window_end] == min(x[peaks[i]:search_window_end], na.rm=T))[1] - 1
        # check amplitude, add segment to list if it covers minimum amplitude. Otherwise, do nothing and move on
        if(x[peaks[i]] - x[potential_trough] > min_seg_amplitude){
          if(is.null(full_segs)){
            full_segs <- list(c(tmp_seg, potential_trough))
          }else{
            full_segs <- c(full_segs, list(c(tmp_seg, potential_trough)))
          }
        }
      }else{
        # so there are troughs, we want to find the one furthest from the peak that yields a segment that is both
        # greater than min_seg_amplitude, and increases the amplitude more than the agg_amp_frac.
        seg_amp <- 0
        true_trough <- NULL
        for(trough in troughs){
          potential_trough <- peaks[i] + trough - 1
          tmp_amp <- x[peaks[i]] - x[potential_trough]
          if(!is.na(tmp_amp)){
            if((tmp_amp >= min_seg_amplitude) & (tmp_amp >= ((1 + agg_amp_frac) * seg_amp))){
              true_trough <- potential_trough
              seg_amp <- tmp_amp
            }
          }
        }
        if(!is.null(true_trough)){
          if(is.null(full_segs)){
            full_segs <- list(c(tmp_seg, true_trough))
          }else{
            full_segs <- c(full_segs, list(c(tmp_seg, true_trough)))
          }
        }
      }# end check for troughs in search window
    }# end check for valid ascending segment
    tmp_seg <- NA # reset tmp_seg to NA
  } # end for
  if(!is.null(full_segs)){
    return(full_segs)
  }else{
    return(NA)
  }
}
#---------------------------------------------------------------------
# returns the index of the first/last value of x that is
# greater/less than the value of thresh. If gup is False (greendown)
# then it returns the first/last value of x that is less/greater than
# the value of thresh. first/last and greater/less determined by first_greater
# NOTE: if thresh is 1 or 0, rounding error can be a problem. Now we round
# the threshold and each of the evi values to 6 decimal places to compensate
GetThresh <- function(thresh_value, x, first_greater=T, gup=T){
  if(gup){
    if(first_greater){
      return(min(which(round(x, 6) >= round(thresh_value, 6))))
    }else{
      return(max(which(round(x, 6) <= round(thresh_value, 6))))
    }
  }else{
    if(first_greater){
      return(min(which(round(x, 6) <= round(thresh_value, 6))))
    }else{
      return(max(which(round(x, 6) >= round(thresh_value, 6))))
    }
  }
}
#---------------------------------------------------------------------
GetSegThresh <- function(seg, x, thresh, gup=T){
  if(gup){
    # check for valid greenup segment
    if(!is.na(seg[1]) & !is.na(seg[2])){
      gup_thresh <- x[seg[1]] + ((x[seg[2]] - x[seg[1]]) * thresh)
      gup_thresh_index <- GetThresh(gup_thresh, x[seg[1]:seg[2]], first_greater=T, gup=T)
      return(gup_thresh_index + seg[1] - 1)
    }else{
      return(NA)
    }
  }else{
    # check for valid greendown segment
    if(!is.na(seg[2]) & !is.na(seg[3])){
      gdown_thresh <- x[seg[3]] + ((x[seg[2]] - x[seg[3]]) * thresh)
      gdown_thresh_index <- GetThresh(gdown_thresh, x[seg[2]:seg[3]], first_greater=T, gup=F)
      return(gdown_thresh_index + seg[2] - 1)
    }else{
      return(NA)
    }
  }
}
#---------------------------------------------------------------------
# example application
# gup_dates <- pred_dates[unlist(lapply(full_segs, GetSegThresh, x, 0.5, T))]
GetPhenoDates <- function(segs, x, dates, gup_threshes, gdown_threshes){
  pheno_dates <- list()
  for(gup_thresh in gup_threshes){
    tmp_dates <- list(dates[unlist(lapply(segs, GetSegThresh, x, gup_thresh, gup=T), use.names=F)])
    if(all(is.na(unlist(tmp_dates, use.names=F)))){
      tmp_dates <- NA
    }
    pheno_dates <- c(pheno_dates, tmp_dates)
  }
  for(gdown_thresh in gdown_threshes){
    tmp_dates <- list(dates[unlist(lapply(segs, GetSegThresh, x, gdown_thresh, gup=F), use.names=F)])
    if(all(is.na(unlist(tmp_dates, use.names=F)))){
      tmp_dates <- NA
    }
    pheno_dates <- c(pheno_dates, tmp_dates)
  }
  return(list(pheno_dates))
}
#---------------------------------------------------------------------
GetSegMetrics <- function(seg, x_smooth, x_raw, smooth_dates, raw_dates, snow_fills){
  if(is.na(seg)){
    return(NA)
  }
  # get the subset of the smoothed and original time series
  tmp_seg_smooth <- x_smooth[seg[1]:seg[3]]
  tmp_gup_smooth <- x_smooth[seg[1]:seg[2]]
  tmp_gdown_smooth <- x_smooth[seg[2]:seg[3]]
  # # get the greenup segment minimum/maximum SVI
  # seg_min <- min(tmp_gup_smooth, na.rm=T)
  # seg_max <- max(tmp_gup_smooth, na.rm=T)
  # get the full segment minimum/maximum SVI
  seg_min_gup <- min(tmp_gup_smooth, na.rm=T)
  seg_max_gup <- max(tmp_gup_smooth, na.rm=T)
  seg_min_gdown <- min(tmp_gdown_smooth, na.rm=T)
  seg_max_gdown <- max(tmp_gdown_smooth, na.rm=T)
  # get the segment integrated SVI: the sum of values above the segment minimum
  # seg_int <- sum(tmp_seg_smooth)
  seg_min <- min(tmp_seg_smooth, na.rm=T)
  seg_max <- max(tmp_seg_smooth, na.rm=T)
  seg_int <- sum(tmp_seg_smooth - seg_min)
  # get greenup segment spline R^2
  gup_raw_date_inds <- which(raw_dates >= smooth_dates[seg[1]] & raw_dates <= smooth_dates[seg[2]]) # indices in raw data of gup segment
  gup_smooth_date_inds <- match(raw_dates[gup_raw_date_inds], smooth_dates) # indices of raw dates in smooth dates
  gup_raw_data <- x_raw[gup_raw_date_inds] # get the raw data associated with the gup segment
  gup_smooth_data <- x_smooth[gup_smooth_date_inds] # get the smoothed values associated with each raw data value
  gup_snow_data <- snow_fills[gup_raw_date_inds] # get the snow data associated with the gup segment
  # calculate the coeff of determination for the spline fit to the greenup segment
  gup_seg_rsquared <- 1 - (sum((gup_raw_data - gup_smooth_data)^2, na.rm=T) / sum((gup_raw_data - mean(gup_raw_data, na.rm=T))^2, na.rm=T))
  # get greenup segment spline R^2
  gdown_raw_date_inds <- which(raw_dates >= smooth_dates[seg[2]] & raw_dates <= smooth_dates[seg[3]]) # indices in raw data of gdown segment
  gdown_smooth_date_inds <- match(raw_dates[gdown_raw_date_inds], smooth_dates) # indices of raw dates in smooth dates
  gdown_raw_data <- x_raw[gdown_raw_date_inds] # get the raw data associated with the gdown segment
  gdown_smooth_data <- x_smooth[gdown_smooth_date_inds] # get the smoothed values associated with each raw data value
  gdown_snow_data <- snow_fills[gdown_raw_date_inds] # get the raw data associated with the gdown segment
  # calculate the coeff of determination for the spline fit to the greenup segment
  gdown_seg_rsquared <- 1 - (sum((gdown_raw_data - gdown_smooth_data)^2, na.rm=T) / sum((gdown_raw_data - mean(gdown_raw_data, na.rm=T))^2, na.rm=T))
  # get the segment missing raw data
  gup_missing <-  mean(diff(c(seg[1],gup_smooth_date_inds,seg[2])))
  gdown_missing <- mean(diff(c(seg[2],gdown_smooth_date_inds,seg[3])))
  # get the segment snow fill information
  gup_longest <- max(diff(c(seg[1],gup_smooth_date_inds,seg[2])))
  gdown_longest <- max(diff(c(seg[2],gdown_smooth_date_inds,seg[3])))
  #gup_snowfilled <- sum(gup_snow_data != 0, na.rm=T)
  #gdown_snowfilled <- sum(gdown_snow_data != 0, na.rm=T)
  # return the metrics as: seg_min, seg_max, seg_int, gup_r2, gup_missing_percentage, gdown_missing_percentage
  # return(c(seg_min, seg_max, seg_int, gup_seg_rsquared, gup_missing, gdown_seg_rsquared, gdown_missing))
  return(c(seg_min_gup, seg_max_gup, seg_min_gdown, seg_max_gdown, seg_int, gup_seg_rsquared, gup_missing, gup_longest, gdown_seg_rsquared, gdown_missing, gdown_longest))
}
#---------------------------------------------------------------------
DoPhenologyLandsat <- function(x, dates, pheno_pars){
  # run the landsat filter/smoother
  filt_evi2 <- x
  snow_fills <- rep(0, length(filt_evi2))
  # LandsatSmoother <- function(x, dates, fill_quant=0.05, x_min=0, spike_thresh=2, min_resid=0.1, spline_spar=NULL, fill_doy=NULL){
  tmp <- LandsatSmoother(x, dates, fill_quant=pheno_pars$LandsatFillQuant, x_min=pheno_pars$LandsatXmin, spike_thresh=pheno_pars$LandsatSpikeThresh, min_resid=pheno_pars$LandsatMinResid, spline_spar=pheno_pars$spline_spar, fill_doy=pheno_pars$LandsatFillDOY, doAnnual=pheno_pars$LandsatDoAnnual, padHeadTail=pheno_pars$LandsatPadHeadTail)
  if(is.na(tmp)){
    # all values were missing
    # NOTE: probably better way to handle this condition than in LandsatSmoother...
    return(NA)
  }
  pred_dates <- tmp[[1]]
  smooth_evi2 <- tmp[[2]]
  # find valid peaks in the time series
  valid_peaks <- try(FilterPeaks(smooth_evi2, FindPotentialPeaks_minmax(smooth_evi2), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance, min_peak_quantile=pheno_pars$min_peak_quantile), silent=T)
  if(inherits(valid_peaks, 'try-error'))
    return(NA)
  # find full segments
  full_segs <- try(GetFullSegs(smooth_evi2, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, agg_amp_frac=pheno_pars$agg_amp_frac), silent=T)
  if(inherits(full_segs, 'try-error')){
    return(NA)
  }
  # get PhenoDates
  pheno_dates <- try(GetPhenoDates(full_segs, smooth_evi2, pred_dates, pheno_pars$gup_threshes, pheno_pars$gdown_threshes), silent=T)
  if(inherits(pheno_dates, "try-error")){
    return(NA)
  }
  # get the segment metrics
  if(!all(is.na(unlist(pheno_dates, use.names=F)))){
    # seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, x, pred_dates, dates)
    # seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, filt_evi2, pred_dates, dates)
    seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, filt_evi2, pred_dates, dates, snow_fills)
    seg_min_gup <- list(unlist(seg_metrics, use.names=F)[seq(1, length(unlist(seg_metrics, use.names=F)), by=11)])
    seg_max_gup <- list(unlist(seg_metrics, use.names=F)[seq(2, length(unlist(seg_metrics, use.names=F)), by=11)])
    seg_min_gdown <- list(unlist(seg_metrics, use.names=F)[seq(3, length(unlist(seg_metrics, use.names=F)), by=11)])
    seg_max_gdown <- list(unlist(seg_metrics, use.names=F)[seq(4, length(unlist(seg_metrics, use.names=F)), by=11)])
    seg_int <- list(unlist(seg_metrics, use.names=F)[seq(5, length(unlist(seg_metrics, use.names=F)), by=11)])
    gup_rsq <- list(unlist(seg_metrics, use.names=F)[seq(6, length(unlist(seg_metrics, use.names=F)), by=11)])
    gup_missing <- list(unlist(seg_metrics, use.names=F)[seq(7, length(unlist(seg_metrics, use.names=F)), by=11)])
    gup_longest <- list(unlist(seg_metrics, use.names=F)[seq(8, length(unlist(seg_metrics, use.names=F)), by=11)])
    #gup_snowfills <- list(unlist(seg_metrics, use.names=F)[seq(6, length(unlist(seg_metrics, use.names=F)), by=9)])
    gdown_rsq <- list(unlist(seg_metrics, use.names=F)[seq(9, length(unlist(seg_metrics, use.names=F)), by=11)])
    gdown_missing <- list(unlist(seg_metrics, use.names=F)[seq(10, length(unlist(seg_metrics, use.names=F)), by=11)])
    gdown_longest <- list(unlist(seg_metrics, use.names=F)[seq(11, length(unlist(seg_metrics, use.names=F)), by=11)])
    #gdown_snowfills <- list(unlist(seg_metrics, use.names=F)[seq(9, length(unlist(seg_metrics, use.names=F)), by=9)])
    # append this new information to the output
    # pheno_dates <- c(unlist(pheno_dates, rec=F), seg_min, seg_max, seg_int, gup_rsq, gup_missing, gdown_rsq, gdown_missing)
    pheno_dates <- c(unlist(pheno_dates, rec=F, use.names=F), seg_min_gup, seg_max_gup, seg_min_gdown, seg_max_gdown, seg_int, gup_rsq, gup_missing, gup_longest, gdown_rsq, gdown_missing, gdown_longest)
  } else{
    pheno_dates <- c(unlist(pheno_dates, rec=F, use.names=F), NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  }
  return(pheno_dates)
  #return(list(pheno_dates,smooth_evi2,pred_dates))
}
#---------------------------------------------------------------------
PlotSeriesLandsat <- function(x, dates, pheno_pars,site,index,ground.obs,pred_dates=NA){
  # PLOT
  x11(h=5,w=11)
  #pdf(h=5,w=11,paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/',site,'_',index,'.pdf',sep=''))
  par(mar=c(3, 4.5, 3, 1))
  # if prediction dates aren't given, assume a daily time series
  # if(is.na(pred_dates)) pred_dates <- min(dates, na.rm=T):max(dates, na.rm=T)
  # get a time series of snow-filtered EVI2
  filt_evi2 <- x
  # filt_evi2 <- try(SnowFilterEVI2(x, pheno_pars$evi2_snow_quant, pheno_pars$ndsi_thresh), silent=T)
  # if(inherits(filt_evi2, 'try-error'))
  # return(NA)
  # tmp <- LandsatScreenSmooth(x, dates)
  tmp <- LandsatSmoother(x, dates, fill_quant=pheno_pars$LandsatFillQuant, spike_thresh=pheno_pars$LandsatSpikeThresh, min_resid=pheno_pars$LandsatMinResid, spline_spar=pheno_pars$spline_spar, fill_doy=pheno_pars$LandsatFillDOY)
  pred_dates <- tmp[[1]]
  smooth_evi2 <- tmp[[2]]
  # PLOT
  plot(dates, filt_evi2, type="n", col=rgb(0.2, 0.2, 0.2), pch=4, cex=1, xlab="", ylab=index, main=site)
  # smooth and eliminate outliers
  # smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_sigma=out_sigma, spline_spar=spline_spar, out_iterations=out_iterations, pred_dates=pred_dates), silent=T)
  # smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)
  # if(inherits(smooth_evi2, 'try-error'))
  # return(NA)
  # PLOT
  # points(pred_dates, smooth_evi2, type="l", col=brewer.pal(5, "Blues")[4], lwd=1.5)
  # find valid peaks in the time series
  valid_peaks <- try(FilterPeaks(smooth_evi2, FindPotentialPeaks_minmax(smooth_evi2), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)
  if(inherits(valid_peaks, 'try-error'))
    return(NA)
  # find full segments
  full_segs <- try(GetFullSegs(smooth_evi2, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, pheno_pars$agg_amp_frac))
  if(inherits(full_segs, 'try-error')){
    return(NA)
  }else if(!is.na(full_segs)){
    # PLOT
    parcoords <- par()$usr #get plot coordinates
    # greenup.color <- "lightgreen"
    # greendown.color <- "pink"
    greenup.color <- brewer.pal(5, "Greens")[3]
    greendown.color <- brewer.pal(5, "Oranges")[3]
    transp=0.5
    greenup.color <- rgb(t(col2rgb(greenup.color))/255, alpha=transp)
    greendown.color <- rgb(t(col2rgb(greendown.color))/255, alpha=transp)
    for(seg in full_segs){
      gup.poly.x <- c(rep(pred_dates[seg[1]], 2), rep(pred_dates[seg[2]], 2))
      gup.poly.y <- c(parcoords[3], parcoords[4], parcoords[4], parcoords[3])
      gdown.poly.x <- c(rep(pred_dates[seg[2]], 2), rep(pred_dates[seg[3]], 2))
      gdown.poly.y <- c(parcoords[3], parcoords[4], parcoords[4], parcoords[3])
      polygon(gup.poly.x, gup.poly.y, col=greenup.color, border=NA)
      polygon(gdown.poly.x, gdown.poly.y, col=greendown.color, border=NA)
    }
  }
  # get PhenoDates
  pheno_dates <- try(GetPhenoDates(full_segs, smooth_evi2, pred_dates, pheno_pars$gup_threshes, pheno_pars$gdown_threshes))
  if(!inherits(pheno_dates, "try-error")){
    coords <- par()$usr
    for(i in 1:length(pheno_pars$gup_threshes)){
      abline(v=pheno_dates[[1]][[i]], lty=2, col="darkgreen")
      text(pheno_dates[[1]][[i]] - 13, y=coords[4] - (0.1 * (coords[4] - coords[3])), labels=as.Date(pheno_dates[[1]][[i]], origin=as.Date("1970-1-1")), srt=90, col="darkgreen", cex=0.7)
    }
    for(i in 1:length(pheno_pars$gdown_threshes)){
      abline(v=pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]], lty=2, col="darkred")
      text(pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]] + 13, y=coords[4] - (0.1 * (coords[4] - coords[3])), labels=as.Date(pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]], origin=as.Date("1970-1-1")), srt=90, col="darkred", cex=0.7)
    }
  }
  points(dates, filt_evi2, type="p", col=rgb(0.2, 0.2, 0.2), pch=16, cex=0.75)
  points(pred_dates, smooth_evi2, type="l", col=brewer.pal(5, "Blues")[4], lwd=1.5)
  
  for (i in 1:nrow(ground.obs)){
    abline(v=ground.obs[i,1],lwd=3,col='darkgreen')
    abline(v=ground.obs[i,2],lwd=3,col='purple')
    abline(v=ground.obs[i,3],lwd=3,col='red') 
  }
  
  #dev.off()
}

#---------------------------------------------------------------------
PlotSeriesLandsat3 <- function(x, dates, pheno_pars, site, index, pred_dates=NA){
  # PLOT
  x11(h=5,w=11)
  #pdf(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/pdfs/',site,'_',index,'.pdf',sep=''))
  par(mar=c(3, 4.5, 3, 1))
  # if prediction dates aren't given, assume a daily time series
  # if(is.na(pred_dates)) pred_dates <- min(dates, na.rm=T):max(dates, na.rm=T)
  # get a time series of snow-filtered EVI2
  filt_evi2 <- x
  # filt_evi2 <- try(SnowFilterEVI2(x, pheno_pars$evi2_snow_quant, pheno_pars$ndsi_thresh), silent=T)
  # if(inherits(filt_evi2, 'try-error'))
  # return(NA)
  # tmp <- LandsatScreenSmooth(x, dates)
  tmp <- LandsatSmoother(x, dates, fill_quant=pheno_pars$LandsatFillQuant, spike_thresh=pheno_pars$LandsatSpikeThresh, min_resid=pheno_pars$LandsatMinResid, spline_spar=pheno_pars$spline_spar, fill_doy=pheno_pars$LandsatFillDOY)
  pred_dates <- tmp[[1]]
  smooth_evi2 <- tmp[[2]]
  # PLOT
  plot(dates, filt_evi2, type="n", col=rgb(0.2, 0.2, 0.2), pch=4, cex=1, xlab="", ylab=index, main=site)
  # smooth and eliminate outliers
  # smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_sigma=out_sigma, spline_spar=spline_spar, out_iterations=out_iterations, pred_dates=pred_dates), silent=T)
  # smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)
  # if(inherits(smooth_evi2, 'try-error'))
  # return(NA)
  # PLOT
  # points(pred_dates, smooth_evi2, type="l", col=brewer.pal(5, "Blues")[4], lwd=1.5)
  # find valid peaks in the time series
  valid_peaks <- try(FilterPeaks(smooth_evi2, FindPotentialPeaks_minmax(smooth_evi2), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)
  if(inherits(valid_peaks, 'try-error'))
    return(NA)
  # find full segments
  full_segs <- try(GetFullSegs(smooth_evi2, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, pheno_pars$agg_amp_frac))
  if(inherits(full_segs, 'try-error')){
    return(NA)
  }else if(!is.na(full_segs)){
    # PLOT
    parcoords <- par()$usr #get plot coordinates
    # greenup.color <- "lightgreen"
    # greendown.color <- "pink"
    greenup.color <- brewer.pal(5, "Greens")[3]
    greendown.color <- brewer.pal(5, "Reds")[3]
    transp=0.5
    greenup.color <- rgb(t(col2rgb(greenup.color))/255, alpha=transp)
    greendown.color <- rgb(t(col2rgb(greendown.color))/255, alpha=transp)
    for(seg in full_segs){
      gup.poly.x <- c(rep(pred_dates[seg[1]], 2), rep(pred_dates[seg[2]], 2))
      gup.poly.y <- c(parcoords[3], parcoords[4], parcoords[4], parcoords[3])
      gdown.poly.x <- c(rep(pred_dates[seg[2]], 2), rep(pred_dates[seg[3]], 2))
      gdown.poly.y <- c(parcoords[3], parcoords[4], parcoords[4], parcoords[3])
      polygon(gup.poly.x, gup.poly.y, col=greenup.color, border=NA)
      polygon(gdown.poly.x, gdown.poly.y, col=greendown.color, border=NA)
    }
  }
  # get PhenoDates
  pheno_dates <- try(GetPhenoDates(full_segs, smooth_evi2, pred_dates, pheno_pars$gup_threshes, pheno_pars$gdown_threshes))
#   if(!inherits(pheno_dates, "try-error")){
#     coords <- par()$usr
#     for(i in 1:length(pheno_pars$gup_threshes)){
#       abline(v=pheno_dates[[1]][[i]], lty=2, col="darkgreen")
#       text(pheno_dates[[1]][[i]] - 13, y=coords[4] - (0.1 * (coords[4] - coords[3])), labels=as.Date(pheno_dates[[1]][[i]], origin=as.Date("1970-1-1")), srt=90, col="darkgreen", cex=0.7)
#     }
#     for(i in 1:length(pheno_pars$gdown_threshes)){
#       abline(v=pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]], lty=2, col="darkred")
#       text(pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]] + 13, y=coords[4] - (0.1 * (coords[4] - coords[3])), labels=as.Date(pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]], origin=as.Date("1970-1-1")), srt=90, col="darkred", cex=0.7)
#     }
#   }
  points(dates, filt_evi2, type="p", col=rgb(0.2, 0.2, 0.2), pch=16, cex=1)
  points(pred_dates, smooth_evi2, type="l", col=brewer.pal(5, "Blues")[4], lwd=2)
  
  #dev.off()

  return(list(smooth_evi2,pred_dates,pheno_dates))
}



#PlotSeriesLandsat(as.matrix(index3[,4]),time,pheno_pars)

#---------------------------------------------------------------------

# ####### Loop through all PhenoCam / Fluxnet locations #######
# setwd('/projectnb/modislc/projects/landsat_sentinel/sites/')
# d <- dir()
# site_info <- rbind(c('HBB','hubbard',43.9438,-71.701,18,5000,'US','nh','18TYP'),
#   c('HVF','harvard',42.5378,-72.1715,18,5000,'US','ma','18TYN'),
#   c('LET','lethbridge',49.709190,-112.940250,11,5000,'CA','NA','11UQR'),
#   c('RUS','russellsage',32.456961,-91.974322,15,5000,'US','la','15SWR'),
#   c('VAI','vaira',38.413281,-120.950636,10,5000,'US','ca','10SFH'))
# 
# site_info <- data.frame(site_info[,1],site_info[,2],
#   as.numeric(site_info[,3]),as.numeric(site_info[,4]),as.numeric(site_info[,5]),as.numeric(site_info[,6]),
#   site_info[,7],site_info[,8],site_info[,9])
# colnames(site_info) <- c('code','pcam','lat','lon','UTM','buffwidth','country','state','granule')
# 
# for (k in 1:5){
#   print(k)
#   site <- site_info[k,]
#   
#   #Find all images for each site
#   setwd('/projectnb/modislc/projects/landsat_sentinel/sites/')
#   
#   VI <- 'EVI1'
#   
#   data_loc <- paste(getwd(),'/',site_info$code[k],'/',VI,'_v1_3',sep='')
#   in_dirs <- list.files(path=data_loc,pattern=glob2rx("HLS**3.tif"),full.names=T,include.dirs=T,recursive=TRUE)
#   tiles <- unique(substr(in_dirs,76,78))
#   in_dirs_tile <- list.files(path=data_loc,pattern=glob2rx(paste("HLS*",tiles[1],"*3.tif",sep="")),
#     full.names=T,include.dirs=T,recursive=TRUE)
#     
#   r <- stack(in_dirs_tile)
#   c <- crop(r,extent(r, 30*(chunk-1)+1, 30*(chunk), 1, 3660))
#   evi_vals <- getValues(c)
#   evi_vals[evi_vals<=0] <- NA
#   
#   #Generate index of ancillary info for all images
#   doy <- as.numeric(substring(in_dirs_tile,85,87))
#   yr <- as.numeric(substring(in_dirs_tile,81,84))
#   sat <- substring(in_dirs_tile,70,70)
#   utm <- as.numeric(substring(in_dirs_tile,75,76))
#   
#   index <- data.frame(sat,yr,doy,utm,t(evi_vals))
#   colnames(index) <- c('sat','yr','doy','utm','evi')
#   index <- index[order(yr,doy),]
#   index$time <- strptime(paste(index$yr, index$doy), format="%Y %j")
#   
#   #index2 <- index2[-which(index2$sat=='S'),]
#   
#   evi_vals <- index[,-c(1,2,3,4,ncol(index))]
#   if (ncol(evi_vals)<10000){
#     block_width <- ncol(evi_vals)-1
#   } else {
#     block_width <- 10000
#   }
#   nblocks <- ncol(evi_vals)%/%block_width
#   bs_start <- seq(1,nblocks*block_width+1,block_width)
#   bs_end <- c(seq(block_width,nblocks*block_width,block_width),ncol(evi_vals))
#   
#   dt.evi <- data.table(evi_vals,keep.rownames=FALSE)
#   info <- as.matrix(index[,2:4])
#   
#   system.time({
#       pheno_mat_all <- foreach(i = 1:ncol(evi_vals), .combine='comb', .multicombine=TRUE,
#         .init=list(list(), list(), list(), list(), list(), list())) %dopar% {
#         if (i%%10000==0) print(i)
#         
#         all <- list(as.numeric(matrix(NA,14,1)),as.numeric(matrix(NA,14,1)),
#           as.numeric(matrix(NA,14,1)),as.numeric(matrix(NA,14,1)),as.numeric(matrix(NA,14,1)))
#           
#         #Generate index table
#         index2 <- cbind(info,dt.evi[,..i]/10000)
#         time <- as.Date(substr(strptime(paste(info[,1], info[,2]), format="%Y %j"),1,10))
#         daynum <- as.numeric(time-time[1])
#         index3 <- index2[-which(is.na(index2[,4])==1),]
#         time <- time[-which(is.na(index2[,4]))]
#           
#         pheno_metrics <- try(DoPhenologyLandsat(as.matrix(index3[,4]),time,pheno_pars))          
# 
#         if (nrow(index3)>0 & is.na(pheno_metrics)==0){
#           num_obs <- nrow(index3) #Calculate number of cloud-free observations
#           w_gup <- which(pheno_metrics[[8]]-pheno_metrics[[7]]>0.1) #which green up shoulders have amplitude > 0.10?
#           w_gdown <- which(pheno_metrics[[10]]-pheno_metrics[[9]]>0.1) #which green down shoulders?
#           
#           num_cyc_gup <- length(w_gup) #total number of green ups
#           num_cyc_gdown <- length(w_gdown) #total number of green downs
#           
#           if (min(c(num_cyc_gup,num_cyc_gdown))>0){
#             all_yrs <- as.numeric(substr(pheno_metrics[[1]],1,4))
#             w_gup_yrs <- all_yrs[w_gup]-2012
#             w_gdown_yrs <- all_yrs[w_gdown]-2012
#             
#             
#             # spr1, spr2, spr3, seg_min_gup, seg_max_gup, gup_rsq, gup_missing, gup_longest
#             tmp <- unlist(pheno_metrics)
#             for (j in 1:length(w_gup)){
#               a <- tmp[seq(w_gup[j],length(tmp),length(all_yrs))] #green up metrics            
#               all[[w_gup_yrs[j]]][c(1:8)] <- a[c(1,2,3,7,8,12,13,14)]
#             }
#             
#             # aut1, aut2, aut3, gdown_rsq, gdown_missing, gdown_longest
#             for (j in 1:length(w_gdown)){
#               b <- tmp[seq(w_gdown[j],length(tmp),length(all_yrs))] #green down metrics
#               all[[w_gdown_yrs[j]]][c(9:14)] <- b[c(4,5,6,15,16,17)]
#             }
#           } else {
#             num_obs <- nrow(index3)
#           }
#           
#         } else {
#           num_obs <- nrow(index3)
#           num_cyc_gup <- 0
#           num_cyc_gdown <- 0
#         }          
#           
#         list(c(num_obs,num_cyc_gup,num_cyc_gdown),all[[1]],all[[2]],all[[3]],all[[4]],all[[5]])
#       }
#       
#   })
#   
#   setwd(paste('/projectnb/modislc/projects/landsat_sentinel/v1_3/Rdata/',site$code,sep=""))
#   save(pheno_mat_all,file=paste(site$code,'_',VI,'_phenology_',chunk,sep=""))
# }

