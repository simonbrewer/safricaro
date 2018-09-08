###############################################################################
## Helper functions for "runSplash.R"
dyn.load('src/splash.so')

###############################################################################
## aetpet: Function to calculate aet and pet
splashf <- function(dtemp,dprec,dsun,lat,yr,elv,whc) {
  dyn.load('src/splash.so')
  retdata = .Fortran('spin_up_sm',
                    yr = as.integer(yr),
                    lat = as.double(lat),
                    elv = as.double(elv),
                    whc = as.double(whc),
                    pr = as.double(dprec),
                    tc = as.double(dtemp),
                    sf = as.double(dsun),
                    maet = double(12),
                    mpet = double(12),
                    mcn = double(12),
                    mro = double(12),
                    mlsr = double(12),
                    msm = double(12),
                    msnow = double(12),
                    daet = double(365),
                    dpet = double(365),
                    dcn = double(365),
                    dro = double(365),
                    dlsr = double(365),
                    dsm = double(365),
                    dsnow = double(365),
                    sm = double(1))
  
  return(retdata)
}
# 
###############################################################################
## aetpet: Function to calculate aet and pet
splashfwbm <- function(dtemp,dprec,dsun,lat,yr,elv,whc) {
  dyn.load('src/splashwbmplus.so')
  retdata = .Fortran('spin_up_sm',
                     yr = as.integer(yr),
                     lat = as.double(lat),
                     elv = as.double(elv),
                     whc = as.double(whc),
                     pr = as.double(dprec),
                     tc = as.double(dtemp),
                     sf = as.double(dsun),
                     maet = double(12),
                     mpet = double(12),
                     mcn = double(12),
                     mro = double(12),
                     mlsr = double(12),
                     msm = double(12),
                     msnow = double(12),
                     mdr = double(12),
                     daet = double(365),
                     dpet = double(365),
                     dcn = double(365),
                     dro = double(365),
                     dlsr = double(365),
                     dsm = double(365),
                     dsnow = double(365),
                     dmelt = double(365),
                     ddr = double(365),
                     sm = double(1))
  
  return(retdata)
}

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
