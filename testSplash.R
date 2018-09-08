## Test splash<->R interface

# Use two adjacent N.American grid boxes
# -90.25, 50.25: SM is ~48121
# -90.25, 50.75: SM is ~29938


require(raster)
dyn.load('fortran/splash.so')
source("helpers.R")

## Climate
## Set climate data files
climdir = "~/Dropbox/CroplandSuitability/CRU climatology/"
tmpfile = "cru_tmp_clim_1961-1990.nc"
prefile = "cru_pre_clim_1961-1990.nc"
cldfile = "cru_cld_clim_1961-1990.nc"

## Get climate data
clim_tmp = stack(paste(climdir,tmpfile,sep=''))
clim_pre = stack(paste(climdir,prefile,sep=''))
clim_cld = stack(paste(climdir,cldfile,sep=''))

tmean = c(0.2,0.9,4.7,9.2,13.5,17.1,19.4,19.2,15.7,10.5,6.0,2.4)
clou =  c(31,38,40,44,49,51,61,60,54,45,31,29)
prec = c(119,112,103,104,112,107,84,83,121,182,176,154)

###############################################################################
## Location 1
## Climate
lon = -90.25
lat = 50.25
yr = 2000
doy = 172
elv = 150

tmean = as.numeric(extract(clim_tmp, cbind(lon,lat)))
clou = as.numeric(extract(clim_cld, cbind(lon,lat)))
prec = as.numeric(extract(clim_pre, cbind(lon,lat)))

dtemp = daily(tmean)$dly
dprec = daily(prec)$dly/(365/12)
dsun = (100-daily(clou)$dly)/100

# retval = .Fortran('evap',
#                   lat = as.double(lat),
#                   doy = as.integer(doy),
#                   elv = as.double(142.0),
#                   yr = as.integer(yr),
#                   sf = as.double(1.0),
#                   tc = as.double(23.0),
#                   sw = as.double(0.9),
#                   cn = double(1),
#                   aet = double(1))

# retval = .Fortran('run_one_day',
#                   lat = as.double(lat),
#                   elv = as.double(142.0),
#                   doy = as.integer(doy),
#                   yr = as.integer(yr),
#                   pr = as.double(5.0),
#                   tc = as.double(23.0),
#                   sf = as.double(1.0),
#                   sm = as.double(75),
#                   ro = double(1))
# #out_evap = evap( lat=37.7, doy=172, elv=142.0, yr=2000, sf=1.0, tc=23.0, sw=0.9 )

# retval = .Fortran('run_one_year',
#                   lat = as.double(lat),
#                   elv = as.double(142.0),
#                   yr = as.integer(yr),
#                   pr = as.double(dprec),
#                   tc = as.double(dtemp),
#                   sf = as.double(dsun),
#                   daet = double(365),
#                   dpet = double(365),
#                   dcn = double(365),
#                   dro = double(365),
#                   dsm = double(365),
#                   sm = as.double(75))

retval1 = .Fortran('spin_up_sm',
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
                   sm = double(1),
                   spin_count = integer(1),
                   diff_sm = double(1))

print(paste("*******************************"))
print(paste("Location 1:"))
print(paste("spin",retval1$spin_count, retval1$diff_sm))
print(paste("sm",sum(retval1$dsm)))
print(paste("aet",sum(retval1$daet)))
print(paste("pet",sum(retval1$dpet)))

###############################################################################
## Location 2
## Climate
lon = -90.25
lat = 50.75
yr = 2000
doy = 172
elv = 150

tmean = as.numeric(extract(clim_tmp, cbind(lon,lat)))
clou = as.numeric(extract(clim_cld, cbind(lon,lat)))
prec = as.numeric(extract(clim_pre, cbind(lon,lat)))

dtemp = daily(tmean)$dly
dprec = daily(prec)$dly/(365/12)
dsun = (100-daily(clou)$dly)/100

# retval = .Fortran('evap',
#                   lat = as.double(lat),
#                   doy = as.integer(doy),
#                   elv = as.double(142.0),
#                   yr = as.integer(yr),
#                   sf = as.double(1.0),
#                   tc = as.double(23.0),
#                   sw = as.double(0.9),
#                   cn = double(1),
#                   aet = double(1))

# retval = .Fortran('run_one_day',
#                   lat = as.double(lat),
#                   elv = as.double(142.0),
#                   doy = as.integer(doy),
#                   yr = as.integer(yr),
#                   pr = as.double(5.0),
#                   tc = as.double(23.0),
#                   sf = as.double(1.0),
#                   sm = as.double(75),
#                   ro = double(1))
# #out_evap = evap( lat=37.7, doy=172, elv=142.0, yr=2000, sf=1.0, tc=23.0, sw=0.9 )

# retval = .Fortran('run_one_year',
#                   lat = as.double(lat),
#                   elv = as.double(142.0),
#                   yr = as.integer(yr),
#                   pr = as.double(dprec),
#                   tc = as.double(dtemp),
#                   sf = as.double(dsun),
#                   daet = double(365),
#                   dpet = double(365),
#                   dcn = double(365),
#                   dro = double(365),
#                   dsm = double(365),
#                   sm = as.double(75))

retval2 = .Fortran('spin_up_sm',
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
                   sm = double(1),
                   spin_count = integer(1),
                   diff_sm = double(1))


print(paste("*******************************"))
print(paste("Location 2:"))
print(paste("spin",retval2$spin_count, retval2$diff_sm))
print(paste("sm",sum(retval2$dsm)))
print(paste("aet",sum(retval2$daet)))
print(paste("pet",sum(retval2$dpet)))
