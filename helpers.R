###############################################################################
## Helper functions for "runSplash.R"
dyn.load('src/splash.so')

###############################################################################
## aetpet: Function to calculate aet and pet
splashf <- function(dtemp,dprec,dsun,lat,yr,elv) {
  dyn.load('src/splash.so')
  retdata <- .Fortran("spin_up_sm",
                      yr = as.integer(yr),
                      lat = as.double(lat),
                      elv = as.double(elv),
                      pr = as.double(dprec),
                      tc = as.double(dtemp),
                      sf = as.double(dsun),
                      maet = double(12),
                      mpet = double(12),
                      mcn = double(12),
                      mro = double(12),
                      msm = double(12),
                      daet = double(365),
                      dpet = double(365),
                      dcn = double(365),
                      dro = double(365),
                      dsm = double(365),
                      sm = double(1))
  return(retdata)
}
# 
###############################################################################
# ## DAILY: Function to interpolate from monthly to daily
# ## Replace with 'approx'?
# 
daily <- function(mly) {
  retdata <- .Fortran("daily",
                      mly = as.double(mly),
                      dly = double(365))
  
  return(retdata)
}
# daily <- function(mly) {
#   ## Time parameters
#   midday = c(16,44,75,105,136,166,197,228,258,289,319,350)
#   daysinmonth = c(31,28,31,30,31,30,31,31,30,31,30,31)
#   hoursinday = seq(1,24)
#   dip = pi/180
#   pid = 180/pi
#   
#   dly = rep(NA, 365)
#   vinc=(mly[1]-mly[12])/31.0
#   dly[350]=mly[12]
#   for (i in 351:365) {
#     dly[i] = dly[(i-1)] + vinc
#   }
#   dly[1]=dly[365]+vinc
#   for (i in 2:15) {
#     dly[i] = dly[(i-1)] + vinc
#   }
#   for (i in 1:11) {
#     vinc=(mly[(i+1)]-mly[i])/(midday[(i+1)]-midday[i])
#     dly[midday[i]]=mly[i]
#     for(j in (midday[i]+1):(midday[i+1]-1)) {
#       dly[j]=dly[j-1]+vinc
#       
#     }
#   }
#   return(dly)
# }
# 

############################################################################
# Matrix manipulation methods
#
# For simplicity we have avoided to create generic functions for 'flip' etc.
# and therefore we have to call the corresponding methods coupled to the
# 'matrix' class explicitly, i.e. flip.matrix().
############################################################################
# Flip matrix (upside-down)
flip.matrix <- function(x) {
  mirror.matrix(rotate180.matrix(x))
}

# Mirror matrix (left-right)
mirror.matrix <- function(x) {
  xx <- as.data.frame(x);
  xx <- rev(xx);
  xx <- as.matrix(xx);
  xx;
}

# Rotate matrix 90 clockworks
rotate90.matrix <- function(x) {
  t(mirror.matrix(x))
}

# Rotate matrix 180 clockworks
rotate180.matrix <- function(x) { 
  xx <- rev(x);
  dim(xx) <- dim(x);
  xx;
}

# Rotate matrix 270 clockworks
rotate270.matrix <- function(x) {
  mirror.matrix(t(x))
}

