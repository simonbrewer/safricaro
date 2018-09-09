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

## Climate
## Set climate data files
climdir = "~/Dropbox/Data/climate/cru_cl_2.00/"
tmpfile = "cru_10min_tmp.nc"
prefile = "cru_10min_pre.nc"
cldfile = "cru_10min_sun.nc"
elvfile = "cru_10min_elv.nc"

## STN grid
stnfile = "~/Dropbox/Data/hydrology/basins/ISLSCP_RIVER_ROUTING/data/river_routing_stn_hdeg/stn_basin_id_hd.asc"
## runoff grid
rofile = "~/Dropbox/Data/hydrology/basins/GRCD/composite/runoff/cmp_ro.grd"

## Output directory
outdir = "./wbm_ro/"
outname = "wbm"

## Get climate data
clim_tmp = stack(paste(climdir,tmpfile,sep=''))
clim_pre = stack(paste(climdir,prefile,sep=''))
clim_cld = stack(paste(climdir,cldfile,sep=''))
clim_elv = raster(paste(climdir,elvfile,sep=''), varname='elv')
clim_elv = clim_elv*1000 ## km -> m

stn_grd = raster(stnfile)
r3 = raster(rofile)

## Crop (test)
myext = extent(c(10,40,-35,-15))
clim_tmp = crop(clim_tmp, myext)
clim_pre = crop(clim_pre, myext)
clim_cld = crop(clim_cld, myext)
stn_grd = crop(stn_grd, myext)
r3 = crop(r3, myext)

clim.crds = SpatialPoints(coordinates(clim_tmp))
ncrds = dim(coordinates(clim.crds))[1]

## Convert climate data to matrices
tmp_mat = extract(clim_tmp, clim.crds)
pre_mat = extract(clim_pre, clim.crds)
cld_mat = extract(clim_cld, clim.crds)
elv_mat = extract(clim_elv, clim.crds)

stn_mat = extract(stn_grd, clim.crds)
stn_mat[which(stn_mat==-88)] <- NA

out.df = list(gdd5 = rep(NA,ncrds),
              tmp = matrix(NA,nrow=ncrds, ncol=12),
              pre = matrix(NA,nrow=ncrds, ncol=12),
              aet = matrix(NA,nrow=ncrds, ncol=12),
              pet = matrix(NA,nrow=ncrds, ncol=12),
              alpha = rep(NA,ncrds),
              cn = matrix(NA,nrow=ncrds, ncol=12),
              sm = matrix(NA,nrow=ncrds, ncol=12),
              runoff = matrix(NA,nrow=ncrds, ncol=12),
              lsr = matrix(NA,nrow=ncrds, ncol=12),
              sw = matrix(NA,nrow=ncrds, ncol=12))

# cl <- makeCluster(4)
# registerDoParallel(cl)
for (i in 1:ncrds) {
  if (i%%100 == 0) {print(paste("Doing:",i,"of",ncrds))}
  
  elv = elv_mat[i]
  # if (!is.na(tmp_mat[i,1]) & !is.na(pre_mat[i,1]) &
  #     !is.na(cld_mat[i,1]) & (elv > -500)){
  if (!is.na(tmp_mat[i,1])) {
    
    ## Set from file
    whc = 150
    # whc = whc_val[whc_mat[i]]
    # if (is.na(whc)) whc = 0
    
    if (whc > 0) {
      dtemp = daily(tmp_mat[i,])$dly
      dprec = daily(pre_mat[i,])$dly/(365/12)
      # dsun = (100-daily(cld_mat[i,])$dly)/100
      dsun = daily(cld_mat[i,])$dly/100
      
      dtemp0 = dtemp-5
      out.df$gdd5[i] = sum(ifelse(dtemp0>0, dtemp0, 0))
      out.df$tmp[i,] = tmp_mat[i,]
      out.df$pre[i,] = pre_mat[i,]
      
      aetpet.df = splashfwbm(dtemp,dprec,dsun,lat=coordinates(clim.crds)[i,2],
                          yr=1975, elv=elv, whc=whc)
      
      out.df$aet[i,] = aetpet.df$maet
      out.df$pet[i,] = aetpet.df$mpet
      out.df$alpha[i] = sum(out.df$aet[i,])/sum(out.df$pet[i,])
      ## ALPHA constraints
      if (out.df$alpha[i]<0) { out.df$alpha[i]=0 }
      if (out.df$alpha[i]>1) { out.df$alpha[i]=1 }
      
      out.df$cn[i,] = aetpet.df$mcn
      ## CN constraints
      out.df$cn[i,] = ifelse(out.df$cn[i,]<0,0,out.df$cn[i,])
      
      out.df$sm[i,] = aetpet.df$msm
      out.df$runoff[i,] = aetpet.df$mro
      out.df$lsr[i,] = aetpet.df$mlsr
      out.df$sw[i,] = aetpet.df$msnow
      
      #if (i == 10991) {stop()}
      
    } else {
      #   print(paste("Missing WHC:", coordinates(clim.crds)[i,1], coordinates(clim.crds)[i,2]))
      out.df$runoff[i,] = out.df$pre[i,]
    }
    #stop()
  }
}
# stopCluster(cl)

## New rasters
r = raster(clim_tmp,1)
gdd5.r = setValues(r, out.df$gdd5)
alpha.r = setValues(r, out.df$alpha)  
#par.r = setValues(r, out.df$par)
for (i in 1:12) {
  if (i == 1) {
    aet.r = setValues(r, out.df$aet[,i], layer=i)  
    pet.r = setValues(r, out.df$pet[,i], layer=i)  
    sw.r = setValues(r, out.df$sw[,i], layer=i)  
    cn.r = setValues(r, out.df$cn[,i], layer=i)  
    runoff.r = setValues(r, out.df$runoff[,i], layer=i)
    lsr.r = setValues(r, out.df$lsr[,i], layer=i)
    
  } else {
    aet.r = stack(aet.r, setValues(r, out.df$aet[,i], layer=i))
    pet.r = stack(pet.r, setValues(r, out.df$pet[,i], layer=i))
    sw.r = stack(sw.r, setValues(r, out.df$sw[,i], layer=i))
    cn.r = stack(cn.r, setValues(r, out.df$cn[,i], layer=i))
    runoff.r = stack(runoff.r, setValues(r, out.df$runoff[,i], layer=i))
    lsr.r = stack(lsr.r, setValues(r, out.df$lsr[,i], layer=i))
    
  }
}

writeRaster(gdd5.r, paste0(outdir,outname,"_gdd.nc"), overwrite=TRUE,
            varname="gdd5",longname="Growing Degree Days",varunit="degree days")
writeRaster(alpha.r, paste0(outdir,outname,"_alpha.nc"), overwrite=TRUE,
            varname="alpha",longname="Moisture Index",varunit="")
writeRaster(aet.r, paste0(outdir,outname,"_aet.nc"), overwrite=TRUE,
            varname="aet",longname="Actual Evapotranspiration",varunit="mm m-1")
writeRaster(pet.r, paste0(outdir,outname,"_pet.nc"), overwrite=TRUE,
            varname="pet",longname="Potential Evapotranspiration",varunit="mm m-1")
writeRaster(sw.r, paste0(outdir,outname,"_sw.nc"), overwrite=TRUE,
            varname="sw",longname="Soil Moisture",varunit="mm m-1")
writeRaster(cn.r, paste0(outdir,outname,"_cn.nc"), overwrite=TRUE,
            varname="cn",longname="Condensation",varunit="mm m-1")
writeRaster(runoff.r, paste0(outdir,outname,"_runoff.nc"), overwrite=TRUE,
            varname="runoff",longname="Runoff",varunit="mm m-1")
writeRaster(lsr.r, paste0(outdir,outname,"_lsr.nc"), overwrite=TRUE,
            varname="lsr",longname="Land Surface Runoff",varunit="mm m-1")

