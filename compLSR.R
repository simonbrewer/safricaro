## Compare runoff values
library(raster)

obs = raster("~/Dropbox/Data/hydrology/basins/GRCD/composite/runoff/cmp_ro.grd")
obs.crds = coordinates(obs)

sm = stack("./sm_ro/sm_lsr.nc")
sm = sum(sm)
sm.crds = coordinates(sm)

wasmod = stack("./wasmod_ro/wasmod_lsr.nc")
wasmod = sum(wasmod)

obs.v = extract(obs, sm.crds)
sm.v = extract(sm, sm.crds)
wasmod.v = extract(wasmod, sm.crds)

plot(obs.v, sm.v)
abline(0,1)
plot(obs.v, wasmod.v)
abline(0,1)
plot(sm.v, wasmod.v)
abline(0,1)
