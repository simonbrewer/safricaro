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

############################################################################
## WASMOD-M model (Widen-Nilsson, 2007, 2009)
wasmodm2 <- function(tmean, prec, pet, aet,
                     Ts = 2, Tm = -2, Ac = 0.5, 
                     Ps = 1e-2, Pf = 1e-3)
{
  
  ## R Script to run the WASMOD-M model
  ## Widen-Nilssen et al., 2007, J. Hydro. 340:105-118
  ## Widen-Nilssen et al., 2009, WRR, 45:W05418
  
  ###############################################################################
  ## Output vectors
  sf = rep(0,12)
  rf = rep(0,12)
  sp = rep(0,12)
  sm = rep(0,12)
  aw = rep(0,12)
  ev = rep(0,12)
  sr = rep(0,12)
  fr = rep(0,12)
  tr = rep(0,12)
  lm = rep(0,12)
  
  ## Initialize snow pack and land moisture
  sp.cur = 5
  lm.cur = 150
  
  ## Parameters for spinup
  lm.sum = 0
  lm.diff = 9999
  diff.eps = 1e-4
  
  mylm = NULL
  while(lm.diff > diff.eps) { ## Outer loop 
    
    for (j in 1:12) {
      ## Snowfall
      sf[j] = prec[j]*(1-(exp((tmean[j] - Ts)/(Ts-Tm)))**2)
      sf[j] = ifelse(sf[j] > 0, sf[j], 0)
      
      ## Rainfall
      rf[j] = prec[j]-sf[j]
      
      ## Melt
      sm[j] = (sp.cur+sf[j])*(1-(exp((Tm - tmean[j])/(Ts-Tm)))**2)
      sm[j] = ifelse(sm[j] > 0, sm[j], 0)
      
      ## Store and update snow pack
      sp[j] = sp.cur + sf[j] - sm[j]
      sp.cur = sp[j]
      
      ## Available water
      aw[j] = lm.cur + rf[j] + sm[j]
      
      ## Evaporation
      # ev[j] = pet[j] * (1 - Ac**(aw[j]/pet[j]))
      # ev[j] = ifelse(ev[j] < aw[j], ev[j], aw[j])
      ev[j] = aet[j]
      ev[j] = ifelse(ev[j] < aw[j], ev[j], aw[j])
      
      ## Slow runoff
      sr[j] = Ps * lm.cur
      
      ## Fast runoff
      fr[j] = Pf * lm.cur * (rf[j] + sm[j])
      
      ## Total runoff
      srfr = sr[j]+fr[j]
      aeev = aw[j]+ev[j]
      tr[j] = ifelse(srfr<aeev,srfr,aeev)
      
      ## Store and update land moisture
      lm[j] = lm.cur + (rf[j] + sm[j] - ev[j] - tr[j])
      lm.cur = lm[j]
    }
    lm.diff = abs(lm.sum - sum(lm))
    lm.sum = sum(lm)
    mylm = c(mylm,lm)
  } ## End spinup
  outl = list(tmean=tmean, 
              prec=prec, 
              pet=pet,
              aet=ev, 
              aw = aw,
              rf=rf, sf=sf, 
              sp=sp, sm=sm,
              sr=sr, fr=fr, 
              tr=tr, lm=lm)
  return(outl)
}

snowmod <- function(Ts = 2, Tm = -2) {
  
  ## Example - figure 7 from Xu
  tt = seq(-10,10,by=0.1)
  pp = rep(5, length(tt))
  sp = rep(7.5, length(tt))
  sf=pp*(1-(exp((tt - Ts)/(Ts-Tm)))**2)
  
  plot(tt, ifelse(sf>=0,sf,0), type='l', ylim=c(0,7.5),
       xlab="Air temperature", ylab="Rate", col=4, lwd=2)
  sm=sp*(1-(exp((Tm - tt)/(Ts-Tm)))**2)
  lines(tt, ifelse(sm>=0,sm,0), col=2, lwd=2)
  legend("right", legend=c("Snow","Melt"), lty=1, lwd=2, col=c(4,2))
  # eqn4a = ((tt-Ts)/(Ts-Tm))
  # eqn4b = ifelse(eqn4a<=0,eqn4a,0)
  # eqn4c = exp(eqn4b)**2
  # eqn4d = pp * (1-eqn4c)
}


