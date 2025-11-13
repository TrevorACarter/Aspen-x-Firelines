#### Insert Beginning of Script here to generate the data to here
## TreeMap
## Cook 2024
## NLCD 2019

#### Projecting and writing
SR <- vect("./Boundary/SouthernRockyBoundary_10kmBuff.shp") ## loading in the study area boundary, may or may not need to re:load based on where in script you are starting

## projecting all of the predictors to make sure they work well together
Aspen_Binary.tif <- project(Aspen_Binary.tif, SR)
DouglasFir_Binary.tif <- project(DouglasFir_Binary.tif, SR)
Grassland_Binary.tif <- project(Grassland_Binary.tif, SR)
Lodgepole_Binary.tif <- project(Lodgepole_Binary.tif, SR)
Oak_Binary.tif <- project(Oak_Binary.tif, SR)
Other_Binary.tif <- project(Other_Binary.tif, SR)
PJ_Binary.tif <- project(PJ_Binary.tif, SR)
Ponderosa_Binary.tif <- project(Ponderosa_Binary.tif, SR)
SF_Binary.tif <- project(SF_Binary.tif, SR)
Shrub_Binary.tif <- project(Shrub_Binary.tif, SR)

Aspen_Binary.tif <- crop(Aspen_Binary.tif, SR)
DouglasFir_Binary.tif <- crop(DouglasFir_Binary.tif, SR)
Grassland_Binary.tif <- crop(Grassland_Binary.tif, SR)
Lodgepole_Binary.tif <- crop(Lodgepole_Binary.tif, SR)
Oak_Binary.tif <- crop(Oak_Binary.tif, SR)
Other_Binary.tif <- crop(Other_Binary.tif, SR)
PJ_Binary.tif <- crop(PJ_Binary.tif, SR)
Ponderosa_Binary.tif <- crop(Ponderosa_Binary.tif, SR)
SF_Binary.tif <- crop(SF_Binary.tif, SR)
Shrub_Binary.tif <- crop(Shrub_Binary.tif, SR)


# may not be necessary
ext(DouglasFir_Binary.tif) <- ext(Aspen_Binary.tif);gc()
ext(Grassland_Binary.tif) <- ext(Aspen_Binary.tif);gc()
ext(Lodgepole_Binary.tif) <- ext(Aspen_Binary.tif);gc()
ext(Oak_Binary.tif) <- ext(Aspen_Binary.tif);gc()
ext(Other_Binary.tif) <- ext(Aspen_Binary.tif);gc()
ext(PJ_Binary.tif) <- ext(Aspen_Binary.tif);gc()
ext(Ponderosa_Binary.tif) <- ext(Aspen_Binary.tif);gc()
ext(SF_Binary.tif) <- ext(Aspen_Binary.tif);gc()
ext(Shrub_Binary.tif) <- ext(Aspen_Binary.tif);gc()

test.stack <- c(Aspen_Binary.tif,DouglasFir_Binary.tif)




