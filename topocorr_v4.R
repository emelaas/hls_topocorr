topocorr_v4 <-function(x, slope, aspect, sunelev, sunazimuth, method="cosine", na.value=NA, GRASS.aspect=FALSE, IL.epsilon=0.000001) {
  # topographic correction for image x based on
  # topography and sun location
  
  # IL.epsilon: if IL == 0, the corrected value is Inf (division by zero)
  # adding a tiny increment eliminates the Inf
  
  ## aspect may be GRASS output: counterclockwise from east
  ## or nonGRASS output: clockwise from north
  ## require the latter for further calculations
  ## because sunazimuth is always measured clockwise from north
  if(GRASS.aspect) {
    aspect <- -1 * aspect + 90
    aspect <- (aspect + 360) %% 360
  }
  
  # all inputs are in degrees, but we need radians
  slope <- (pi/180) * slope
  aspect <- (pi/180) * aspect
  sunzenith <- (pi/180) * (90 - sunelev)
  sunazimuth <- (pi/180) * sunazimuth
  
  x.orig <- x
  x[x == na.value] <- NA
  
  IL <- cos(slope) * cos(sunzenith) + sin(slope) * sin(sunzenith) * cos(sunazimuth - aspect)
  IL[IL == 0] <- IL.epsilon
  
  METHODS <- c("none","ccorrection", "scsc", "empirical", "rotational", "illumination")
  method <- pmatch(method, METHODS)
  if (is.na(method)) 
    stop("invalid method")
  if (method == -1) 
    stop("ambiguous method")
  
  if(method == 1){
    ## No correction
    xout <- x
  }
  else if(method == 2) {
    ## C correction
    ## Teillet, Guindon, and Goodenough 1982
    if (length(x)>1000) {r <- sample(length(x),1000)} else {r <- length(x)}
    band.lm <- lm(as.vector(x[r]) ~ as.vector(IL[r]))
    C <- coefficients(band.lm)[[1]]/coefficients(band.lm)[[2]]
    
    xout <- x * (cos(sunzenith) + C) / (IL + C)
  }
  else if(method == 3) {
    ## SCS + S
    ## Soenen et al. 2005
    if (length(x)>1000) {r <- sample(length(x),1000)} else {r <- length(x)}
    band.lm <- lm(as.vector(x[r]) ~ as.vector(IL[r]))
    C <- coefficients(band.lm)[[1]]/coefficients(band.lm)[[2]]
    
    xout <- x * (cos(slope)*cos(sunzenith) + C) / (IL + C)
  }
  else if(method == 4) {
    ## Empirical Statistical (see Hantson & Chuvieco, 2011 IJAEOG)
    ## Teillet et al. 1982
    if (length(x)>1000) {r <- sample(length(x),1000)} else {r <- length(x)}
    band.lm <- lm(as.vector(x[r]) ~ as.vector(IL[r]))
    B <- coefficients(band.lm)[[1]]
    A <- coefficients(band.lm)[[2]]
    xmean <- mean(as.vector(x), na.rm=TRUE)
    
    xout <- x - (A * IL + B) + xmean
  }
  else if(method == 5) {
    ## Rotation Model
    ## Tan et al. 2010
    if (length(x)>1000) {r <- sample(length(x),1000)} else {r <- length(x)}
    band.lm <- lm(as.vector(x[r]) ~ as.vector(IL[r]))
    A <- coefficients(band.lm)[[2]]
    B <- coefficients(band.lm)[[1]]
    
    xout <- x - A * (IL - cos(sunzenith))
  }
  
  else if(method == 6) {
    ## illumination only
    xout <- IL
  }
  
  ## if slope is zero, reflectance does not change
  if(method != 6) 
    xout[slope == 0 & !is.na(slope)] <- x[slope == 0 & !is.na(slope)]
  
  ## if x was a SpatialGridDataFrame, return an object of the same class
  if(class(x.orig) == "SpatialGridDataFrame") {
    x.orig@data[,1] <- as.vector(xout)
    xout <- x.orig
  }
  
  xout
}