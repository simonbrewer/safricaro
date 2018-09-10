## Compare PET values
require(raster)

obs = read.fwf("./Data/Eo150.clim", 
               widths = rep(8,14), col.names = c("lon","lat",month.abb))
head(obs)
coordinates(obs) <- ~lon+lat

r = raster(res=0.5)

pet.obs = rasterize(obs, r, field=month.abb)
plot(pet.obs)
writeRaster(pet.obs, "./Data/pet.obs.mon.nc", overwrite=TRUE)

pet.obs.ann = sum(pet.obs)
plot(pet.obs.ann)
writeRaster(pet.obs.ann, "./Data/pet.obs.ann.nc", overwrite=TRUE)
