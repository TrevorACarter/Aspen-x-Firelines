#### Begin - Loading Dependencies, Set Working Directory ####
library(terra) ## needed for all spatial analysis
setwd("D:/Aspen Firelines/data") ## set working directory to the data folder

#### LandCover ####
SR <- vect("./Boundary/SouthernRockyBoundary_10kmBuff.shp") ## loading in the study area boundary, may or may not need to re:load based on where in script you are starting
Aspen <- rast("./Preprocessed Rasters/Land Cover/s2aspen_prob_10m_binOpt.tif")
NLCD <- rast("./Preprocessed Rasters/Land Cover/NLCD_2019_SR.tif")
TreeMap <- rast("./Preprocessed Rasters/Land Cover/TreeMap2016.tif")

### Projecting and cropping
Aspen <- project(Aspen, SR)
Aspen <- crop(Aspen, SR)

NLCD <- project(NLCD, SR)
NLCD <- crop(NLCD, SR)

TreeMap <- crop(TreeMap, SR)
TreeMap <- project(TreeMap, SR)

ext(Aspen) <- ext(SR);gc()
ext(NLCD) <- ext(SR);gc()
ext(TreeMap) <- ext(SR);gc()

## need to resample the TreeMap (NLCD and Aspen are already 10x10 m)
TreeMap_resample <- resample(TreeMap,Aspen, method = "near")

plot(Aspen);plot(SR, add = TRUE)
plot(NLCD);plot(SR, add = TRUE)
plot(TreeMap_resample);plot(SR, add = TRUE)

test.stack <- c(Aspen,NLCD,TreeMap_resample)
rm(test.stack)

## making the binary rasters for presence/absence
Grass_Binary <- ifel(NLCD == 71,1,0) ## NLCD value of 71 == Herbaceous/Graminoid cover
plot(Grass_Binary)
Shrub_Binary <- ifel(NLCD == 52,1,0) ## NLCD value of 71 == Shrubland cover
plot(Shrub_Binary)

TreeMap_key <- read.csv("./Preprocessed Rasters/Land Cover/TreeMap2016_tree_table.csv")
table(TreeMap_key$COMMON_NAME)

Pinyon <- unique(TreeMap_key$tm_id[grepl("pinyon", TreeMap_key$COMMON_NAME)])
Juniper <- unique(TreeMap_key$tm_id[grepl("juniper", TreeMap_key$COMMON_NAME)])
PJ_Binary <- ifel(TreeMap_resample %in% Pinyon | TreeMap_resample %in% Juniper, 1, 0)
plot(PJ_Binary)
rm(Pinyon);rm(Juniper)

Spruce <- unique(TreeMap_key$tm_id[TreeMap_key$COMMON_NAME == "Engelmann spruce"])
Fir <- unique(TreeMap_key$tm_id[TreeMap_key$COMMON_NAME == "subalpine fir"])
SF_Binary <- ifel(TreeMap_resample %in% Spruce | TreeMap_resample %in% Fir, 1, 0)
plot(SF_Binary)
rm(Spruce);rm(Fir)

Lodgepole <- unique(TreeMap_key$tm_id[TreeMap_key$COMMON_NAME == "lodgepole pine"])
Lodgepole_Binary <- ifel(TreeMap_resample %in% Lodgepole, 1, 0)
plot(Lodgepole_Binary)
rm(Lodgepole)

Ponderosa <- unique(TreeMap_key$tm_id[TreeMap_key$COMMON_NAME == "ponderosa pine"])
Ponderosa_Binary <- ifel(TreeMap_resample %in% Ponderosa, 1, 0)
plot(Ponderosa_Binary)
rm(Ponderosa)

Gambel <- unique(TreeMap_key$tm_id[TreeMap_key$COMMON_NAME == "Gambel oak"])
Gambel_Binary <- ifel(TreeMap_resample %in% Gambel, 1, 0)
plot(Gambel_Binary)
rm(Gambel)

DougFir <- unique(TreeMap_key$tm_id[TreeMap_key$COMMON_NAME == "Douglas-fir"])
DougFir_Binary <- ifel(TreeMap_resample %in% DougFir, 1, 0)
plot(DougFir_Binary)
rm(DougFir)

## other layer (described in the supplemental materials)
stack <- c(Aspen,PJ_Binary,SF_Binary,Lodgepole_Binary,Shrub_Binary,Grass_Binary,Ponderosa_Binary,Gambel_Binary,DougFir_Binary)
all_equal_0 <- app(stack, fun=function(x) all(x == 0))
plot(all_equal_0)

#### writing rasters
writeRaster(Aspen,"./Predictor Rasters/Land Cover/Aspen_Binary.tif", overwrite = TRUE)
writeRaster(PJ_Binary,"./Predictor Rasters/Land Cover/PJ_Binary.tif", overwrite = TRUE)
writeRaster(SF_Binary,"./Predictor Rasters/Land Cover/SF_Binary.tif", overwrite = TRUE)
writeRaster(Lodgepole_Binary,"./Predictor Rasters/Land Cover/Lodgepole_Binary.tif", overwrite = TRUE)
writeRaster(Shrub_Binary,"./Predictor Rasters/Land Cover/Shrub_Binary.tif", overwrite = TRUE)
writeRaster(Grass_Binary,"./Predictor Rasters/Land Cover/Grass_Binary.tif", overwrite = TRUE)
writeRaster(Ponderosa_Binary,"./Predictor Rasters/Land Cover/Ponderosa_Binary.tif", overwrite = TRUE)
writeRaster(Gambel_Binary,"./Predictor Rasters/Land Cover/Gambel_Binary.tif", overwrite = TRUE)
writeRaster(DougFir_Binary,"./Predictor Rasters/Land Cover/DougFir_Binary.tif", overwrite = TRUE)
writeRaster(all_equal_0,"./Predictor Rasters/Land Cover/Other_Binary.tif", overwrite = TRUE)


#### Topographic ####
SR <- vect("./Boundary/SouthernRockyBoundary_10kmBuff.shp") ## loading in the study area boundary, may or may not need to re:load based on where in script you are starting
elev <- rast("./Preprocessed Rasters/Topographic/SRockies_DEM_10kmBuff.tif")
slope <- rast("./Preprocessed Rasters/Topographic/SRockies_slope_10kmBuff.tif")
tpi <- rast("./Preprocessed Rasters/Topographic/SRockies_TPI_10kmBuff.tif")

## projecting and cropping to the study area
elev <- project(elev, SR)
elev <- crop(elev, SR)
ext(elev) <- ext(SR);gc()

slope <- project(slope, SR)
slope <- crop(slope, SR)
ext(slope) <- ext(SR);gc()

tpi <- project(tpi, SR)
tpi <- crop(tpi, SR)
ext(tpi) <- ext(SR);gc()

elev ## looking at the raster to see if it needs additional correction
slope
tpi
aspen <- rast("./Predictor Rasters/Land Cover/Aspen_Binary.tif") 

plot(elev);plot(SR, add = TRUE)
plot(slope);plot(SR, add = TRUE)
plot(tpi);plot(SR, add = TRUE)

## need to resample to match aspen layer directly
elev <- resample(elev,aspen, method = "bilinear")
slope <- resample(slope,aspen, method = "bilinear")
tpi <- resample(tpi,aspen, method = "bilinear")

test.stack <- c(elev,slope,tpi,aspen)
rm(test.stack)

## writing updated rasters
writeRaster(elev,"./Predictor Rasters/Topographic/elev.tif", overwrite = TRUE)
writeRaster(slope,"./Predictor Rasters/Topographic/slope.tif", overwrite = TRUE)
writeRaster(tpi,"./Predictor Rasters/Topographic/tpi.tif", overwrite = TRUE)


#### Climate ####
aspen <- rast("./Predictor Rasters/Land Cover/Aspen_Binary.tif") ## loading in aspen for future resampling
SR <- vect("./Boundary/SouthernRockyBoundary_10kmBuff.shp") ## loading in the study area boundary, may or may not need to re:load based on where in script you are starting
all_files <- list.files(path = "./Preprocessed Rasters/Climate/",pattern="*.tif") ## creating a vector that has all the files in the working directory with .tif extensions

## if computing environment is able to process large datasets, use this line
# temp <- all_files


## if computing environment is not able to process large datasets, split into subsets and save intermediate files
avgRasts <- all_files[grepl("mean", all_files)]
yr2019 <- all_files[grepl("2019", all_files)]
# yr2020 <- all_files[grepl("2020", all_files)]
# yr2021 <- all_files[grepl("2021", all_files)]
# yr2022 <- all_files[grepl("2022", all_files)]
# yr2023 <- all_files[grepl("2023", all_files)]
temp <- yr2019 ## replace will which subset of interest


## for either all data, or subsets, load in the data
for(i in 1:length(temp)) {
  path <- paste("./Preprocessed Rasters/Climate/", temp[i], sep = "") ## specifying the relative pathway for assign
  assign(temp[i], terra::rast(path)) ## assigning the raster
} ## loading in the shapefiles I want

temp

# reprojecting and cropping raster files
raster_list <- mget(temp, envir = .GlobalEnv)
reprojected_list <- lapply(raster_list, function(x) {
  reprojected <- project(x, crs(SR)) ## projecting to the SR crs
  matched <- crop(reprojected, SR) ## cropping to the SR
  ext(matched) <- ext(SR) ## matching the extent
  resampled <- resample(matched,aspen, method = "bilinear") ## resampling to the aspen layer (takes a while)
  return(resampled)
})
list2env(reprojected_list, envir = .GlobalEnv)
gc()


## if choosing subset, write the rasters into an intermediate folder after projecting
output_path <- "./Preprocessed Rasters/Climate_Intermediate/"
lapply(temp, function(name) {
  rast_obj <- get(name, envir = .GlobalEnv)
  output_file <- file.path(output_path, name)
  writeRaster(rast_obj, output_file, overwrite = TRUE)
})

## reading in the intermediate files post projection (memory no longer an issue after projection)
temp <- list.files(path = "./Preprocessed Rasters/Climate_Intermediate/",pattern="*.tif") ## creating a vector that has all the files in the working directory with .tif extensions
for(i in 1:length(temp)) {
  path <- paste("./Preprocessed Rasters/Climate_Intermediate/", temp[i], sep = "") ## specifying the relative pathway for assign
  assign(temp[i], terra::rast(path)) ## assigning the raster
} ## loading in the shapefiles I want

temp

## correcting the units for all of the annual data 
annual_rasters <- temp[grepl("20", temp)]
raster_list <- mget(annual_rasters, envir = .GlobalEnv)
divided_list <- lapply(raster_list, function(x) {
  x /10  # or whatever value you need
})
list2env(divided_list, envir = .GlobalEnv)


## Creating the difference of mean - month
SearchNames <- c("03_20","04_20","05_20","06_20","07_20","08_20","09_20","10_20","11_20") ## all of the annual data share this pattern (by month)
MeanNames <- c("03_mean","04_mean","05_mean","06_mean","07_mean","08_mean","09_mean","10_mean","11_mean")

for(i in 1:length(SearchNames)){
Search_rasters <- temp[grepl(SearchNames[i], temp)]
raster_list <- mget(Search_rasters, envir = .GlobalEnv)
subtract_list <- lapply(raster_list, function(x) {
  Sub_raster <- mget(temp[grepl(MeanNames[i]), temp], envir = .GlobalEnv)
  x - Sub_raster  # or whatever month
})
list2env(subtract_list, envir = .GlobalEnv)
}
gc()
## clean the global environment
## can only have the rasters that are getting written in GlobalEnv?

## write the rasters
output_path <- "./Predictor Rasters/Climate/"
lapply(temp, function(name) {
  rast_obj <- get(name, envir = .GlobalEnv)
  output_file <- file.path(output_path, paste0(name, ".tif"))
  writeRaster(rast_obj, output_file, overwrite = TRUE)
})