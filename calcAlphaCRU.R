require(raster)
require(ncdf4)
source("helpers.R")

## Climate
## Set climate data files
climdir = "~/Dropbox/CroplandSuitability/Climatology/"
tmpfile = "cru_tmp_clim_1961-1990.nc"
prefile = "cru_pre_clim_1961-1990.nc"
cldfile = "cru_cld_clim_1961-1990.nc"

## Get climate data
clim_tmp = stack(paste(climdir,tmpfile,sep=''))
clim_pre = stack(paste(climdir,prefile,sep=''))
clim_cld = stack(paste(climdir,cldfile,sep=''))

## Crop
# myext = extent(c(10,12,64,66))
# clim_tmp = crop(clim_tmp, myext)
# clim_pre = crop(clim_pre, myext)
# clim_cld = crop(clim_cld, myext)

clim.crds = coordinates(clim_tmp)
ncrds = dim(clim.crds)[1]

## Convert climate data to matrices
tmp_mat = getValues(clim_tmp)
pre_mat = getValues(clim_pre)
cld_mat = getValues(clim_cld)

out.df = data.frame(gdd5 = rep(NA,ncrds),
                    aet = rep(NA,ncrds),
                    pet = rep(NA,ncrds),
                    alpha = rep(NA,ncrds),
                    sm = rep(NA,ncrds),
                    runoff = rep(NA,ncrds),
                    snowpack = rep(NA,ncrds))

for (i in 1:ncrds) {
  if (i%%100 == 0) {print(paste("Doing:",i,"of",ncrds))}
  
  if (!is.na(tmp_mat[i,1]) & !is.na(pre_mat[i,1]) & !is.na(cld_mat[i,1])){
    
    dtemp = daily(tmp_mat[i,])$dly
    dprec = daily(pre_mat[i,])$dly/(365/12)
    dsun = (100-daily(cld_mat[i,])$dly)/100
    
    dtemp0 = dtemp-5
    out.df$gdd5[i] = sum(ifelse(dtemp0>0, dtemp0, 0))
    
    aetpet.df = splashf(dtemp,dprec,dsun,lat=clim.crds[i,2],
                        yr=2000, elv=150)
    out.df$aet[i] = sum(aetpet.df$daet, na.rm = TRUE)
    out.df$pet[i] = sum(aetpet.df$dpet, na.rm = TRUE)
    out.df$alpha[i] = out.df$aet[i]/out.df$pet[i]
    out.df$sm[i] = sum(aetpet.df$dsm, na.rm = TRUE)
    out.df$runoff[i] = sum(aetpet.df$dro, na.rm = TRUE)
    #out.df$snowpack[i] = max(aetpet.df$snowpack, na.rm = TRUE)
    #stop()
  }
}

## New rasters
r = raster(clim_tmp,1)
gdd5.r = setValues(r, out.df$gdd5)
#par.r = setValues(r, out.df$par)
aet.r = setValues(r, out.df$aet)
pet.r = setValues(r, out.df$pet)
alpha.r = setValues(r, out.df$alpha)
sw.r = setValues(r, out.df$sm)
runoff.r = setValues(r, out.df$runoff)

## Output
## Get output parameters from existing file
ncid = nc_open(paste(climdir,tmpfile,sep=''))
mylon_att = ncatt_get(ncid, "longitude")
mylon = ncvar_get(ncid, "longitude")
mylat_att = ncatt_get(ncid, "latitude")
mylat = ncvar_get(ncid, "latitude")
nc_close(ncid)

## Dimension definition
dimX <- ncdim_def( "lon", mylon_att$units, mylon, 
                   create_dimvar=TRUE, longname="longitude" )
dimY <- ncdim_def( "lat", mylat_att$units, mylat, 
                   create_dimvar=TRUE, longname="latitude" )

## Variable definition
gdd5_def = ncvar_def("gdd5", units="degdays", list(dimX,dimY), 1.e30 )
pet_def = ncvar_def("pet", units="mm", list(dimX,dimY), 1.e30 )
aet_def = ncvar_def("aet", units="mm", list(dimX,dimY), 1.e30 )
alpha_def = ncvar_def("alpha", units="ratio", list(dimX,dimY), 1.e30 )
sm_def = ncvar_def("sm", units="mm", list(dimX,dimY), 1.e30 )
runoff_def = ncvar_def("runoff", units="mm", list(dimX,dimY), 1.e30 )
snowpack_def = ncvar_def("snowpack", units="mm", list(dimX,dimY), 1.e30 )

ncnew <- nc_create( "bioclim_splash.nc", list(gdd5_def, aet_def, pet_def,
                                       alpha_def, sm_def, runoff_def,
                                       snowpack_def))
ncvar_put(ncnew, gdd5_def, 
          mirror.matrix(matrix(out.df$gdd5, nrow=length(mylon))))
ncvar_put(ncnew, aet_def, 
          mirror.matrix(matrix(out.df$aet, nrow=length(mylon))))
ncvar_put(ncnew, pet_def, 
          mirror.matrix(matrix(out.df$pet, nrow=length(mylon))))
ncvar_put(ncnew, alpha_def, 
          mirror.matrix(matrix(out.df$alpha, nrow=length(mylon))))
ncvar_put(ncnew, sm_def, 
          mirror.matrix(matrix(out.df$sm, nrow=length(mylon))))
ncvar_put(ncnew, runoff_def, 
          mirror.matrix(matrix(out.df$runoff, nrow=length(mylon))))
ncvar_put(ncnew, snowpack_def, 
          mirror.matrix(matrix(out.df$snowpack*10, nrow=length(mylon))))

ncatt_put(ncnew,"gdd5","standard_name","Growing Degree Days",prec="text")
ncatt_put(ncnew,"aet","standard_name","Actual Evapotranspiration",prec="text")
ncatt_put(ncnew,"pet","standard_name","Potential Evapotranspiration",prec="text")
ncatt_put(ncnew,"alpha","standard_name","Moisture Index",prec="text")
ncatt_put(ncnew,"sm","standard_name","Soil Moisture",prec="text")
ncatt_put(ncnew,"runoff","standard_name","Runoff",prec="text")
ncatt_put(ncnew,"snowpack","standard_name","Snowpack",prec="text")

ncatt_put(ncnew,"lon","standard_name","longitude",prec="text")
ncatt_put(ncnew,"lat","standard_name","latitude",prec="text")

nc_close(ncnew)

