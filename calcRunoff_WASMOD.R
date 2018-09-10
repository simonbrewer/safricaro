###############################################################################
## Version 0.6: 
##
## Corrected snow melt in fortran code
## Add elevation
require(foreach)
require(doParallel)
require(raster)
require(ncdf4)
source("helpers.R")

## Parameter set for WASMOD-M
Ts = 2 ## Upper temperature for snow
Tm = -2 ## Lower temperature for snow
Ac = 0.5 ## Actual evapotranspiration
Ps = 1e-3 ## Slow runoff
Pf = 1e-2 ## Fast runoff
ec = 0.018 ## (mm m-1 C-2)

## Climate
## Set climate data files
climdir = "~/Dropbox/Data/climate/cru_cl_2.00/"
tmpfile = "cru_10min_tmp.nc"
prefile = "cru_10min_pre.nc"
cldfile = "cru_10min_sun.nc"
rehfile = "cru_10min_reh.nc"
elvfile = "cru_10min_elv.nc"

## STN grid
stnfile = "~/Dropbox/Data/hydrology/basins/ISLSCP_RIVER_ROUTING/data/river_routing_stn_hdeg/stn_basin_id_hd.asc"
## runoff grid
rofile = "~/Dropbox/Data/hydrology/basins/GRCD/composite/runoff/cmp_ro.grd"

## Output directory
outdir = "./wasmod_ro/"
outname = "wasmod"

## Get climate data
clim_tmp = stack(paste(climdir,tmpfile,sep=''))
clim_pre = stack(paste(climdir,prefile,sep=''))
clim_cld = stack(paste(climdir,cldfile,sep=''))
clim_reh = stack(paste(climdir,rehfile,sep=''))
clim_elv = raster(paste(climdir,elvfile,sep=''), varname='elv')
clim_elv = clim_elv*1000 ## km -> m

stn_grd = raster(stnfile)
r3 = raster(rofile)

## Crop (test)
myext = extent(c(10,40,-35,-15))
clim_tmp = crop(clim_tmp, myext)
clim_pre = crop(clim_pre, myext)
clim_cld = crop(clim_cld, myext)
clim_reh = crop(clim_cld, myext)
stn_grd = crop(stn_grd, myext)
r3 = crop(r3, myext)

clim.crds = SpatialPoints(coordinates(clim_tmp))
ncrds = dim(coordinates(clim.crds))[1]

## Convert climate data to matrices
tmp_mat = extract(clim_tmp, clim.crds)
pre_mat = extract(clim_pre, clim.crds)
cld_mat = extract(clim_cld, clim.crds)
reh_mat = extract(clim_reh, clim.crds)
elv_mat = extract(clim_elv, clim.crds)

stn_mat = extract(stn_grd, clim.crds)
stn_mat[which(stn_mat==-88)] <- NA

out.df = list(gdd5 = rep(NA,ncrds),
              tmp = matrix(NA,nrow=ncrds, ncol=12),
              pre = matrix(NA,nrow=ncrds, ncol=12),
              aet = matrix(NA,nrow=ncrds, ncol=12),
              pet = matrix(NA,nrow=ncrds, ncol=12),
              ep = matrix(NA,nrow=ncrds, ncol=12),
              alpha = rep(NA,ncrds),
              cn = matrix(NA,nrow=ncrds, ncol=12),
              sm = matrix(NA,nrow=ncrds, ncol=12),
              runoff = matrix(NA,nrow=ncrds, ncol=12),
              lsr = matrix(NA,nrow=ncrds, ncol=12),
              sw = matrix(NA,nrow=ncrds, ncol=12))

for (i in 1:ncrds) {
  if (i%%100 == 0) {print(paste("Doing:",i,"of",ncrds))}
  
  elv = elv_mat[i]
  # if (!is.na(tmp_mat[i,1]) & !is.na(pre_mat[i,1]) &
  #     !is.na(cld_mat[i,1]) & (elv > -500)){
  if (!is.na(tmp_mat[i,1])) {
    
    ## Set from file
    whc = 150
    
    if (whc > 0) {

      dtemp0 = dtemp-5
      out.df$gdd5[i] = sum(ifelse(dtemp0>0, dtemp0, 0))
      out.df$tmp[i,] = tmp_mat[i,]
      out.df$pre[i,] = pre_mat[i,]
      out.df$ep[i,] = ec * ifelse(tmp_mat[i,]>0, tmp_mat[i,]^2, 0) * (100 - reh_mat[i,])

      ###############################################################################
      ## Uses original vectors - can be cleaned
      tmean = out.df$tmp[i,]
      prec = out.df$pre[i,]
      pet = out.df$ep[i,]
      
      ## Temporary vectors
      sp.cur = 5
      lm.cur = 150
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
      
      ## Convergence loop
      lm.sum = 0
      lm.diff = 9999
      diff.eps = 1e-4
      i=0
      mylm = NULL
      
      while(lm.diff > diff.eps) { ## Outer loop 
        # print(i)
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
          ev[j] = pet[j] * (1 - Ac**(aw[j]/pet[j]))
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
        } ## Monthly loop
        lm.diff = abs(lm.sum - sum(lm))
        lm.sum = sum(lm)
        mylm = c(mylm,lm)
        i = i+1
      } ## Convergence loop
      out.df$lsr[i,] = tr
      out.df$sm[i,] = sm
      out.df$sw[i,] = sf
    }
  }
  
}

## New rasters
r = raster(clim_tmp,1)
gdd5.r = setValues(r, out.df$gdd5)
#par.r = setValues(r, out.df$par)
for (i in 1:12) {
  if (i == 1) {
    pet.r = setValues(r, out.df$ep[,i], layer=i)  
    sm.r = setValues(r, out.df$sm[,i], layer=i)  
    sw.r = setValues(r, out.df$sw[,i], layer=i)  
    lsr.r = setValues(r, out.df$lsr[,i], layer=i)
    
  } else {
    pet.r = stack(pet.r, setValues(r, out.df$ep[,i], layer=i))
    sm.r = stack(sm.r, setValues(r, out.df$sm[,i], layer=i))
    sw.r = stack(sw.r, setValues(r, out.df$sw[,i], layer=i))
    lsr.r = stack(lsr.r, setValues(r, out.df$lsr[,i], layer=i))
    
  }
}

writeRaster(gdd5.r, paste0(outdir,outname,"_gdd.nc"), overwrite=TRUE,
            varname="gdd5",longname="Growing Degree Days",varunit="degree days")
writeRaster(pet.r, paste0(outdir,outname,"_pet.nc"), overwrite=TRUE,
            varname="pet",longname="Potential Evapotranspiration",varunit="mm m-1")
writeRaster(sm.r, paste0(outdir,outname,"_sm.nc"), overwrite=TRUE,
            varname="sm",longname="Soil moisture",varunit="mm m-1")
writeRaster(sw.r, paste0(outdir,outname,"_sw.nc"), overwrite=TRUE,
            varname="sw",longname="Snow pack",varunit="mm")
writeRaster(lsr.r, paste0(outdir,outname,"_lsr.nc"), overwrite=TRUE,
            varname="lsr",longname="Land Surface Runoff",varunit="mm m-1")

