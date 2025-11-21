#### Begin - Loading Dependencies, Set Working Directory ####
library(terra) ## needed for all spatial analysis
library(exactextractr) ## faster extract function than terra::extract
library(dplyr) ## less error associated with rounding compared to unique()
library(ncf) ## needed to compute correlograms for spatial autocorrelation
library(vegan) ## needed for the PCNM method to account for spatial autocorrelation
library(RSpectra) ## faster way to calculate eigen vectors
library(caret) ## for random forest training
library(randomForest) ## for random forest models
library(pROC) ## needed for AUC
library(randomForestSRC) ## for tuning RF
setwd("D:/Aspen Firelines/data") ## set working directory to the data folder


#### Processing Fire Lines Data from NIFC (2019-2023) ####
temp <- list.files(path = "./NIFC Lines/",pattern="*.shp") ## creating a vector that has all the files in the working directory with .shp extensions
temp <- temp[1:5] ## setting 1:5 as script will write a compiled .shp file (stored in archived data)

for(i in 1:length(temp)) {
  path <- paste("./NIFC Lines/", temp[i], sep = "") ## specifying the relative pathway for assign
  assign(temp[i], terra::vect(path)) ## assigning the shapefiles
  } ## loading in the shapefiles I want

EventLine2019.shp$year <- 2019
EventLine2020.shp$year <- 2020
EventLine2021.shp$year <- 2021
EventLine2022.shp$year <- 2022
EventLine2023.shp$year <- 2023

stacked_FL <- rbind(EventLine2019.shp,EventLine2020.shp,EventLine2021.shp,EventLine2022.shp,EventLine2023.shp)
rm(EventLine2019.shp);rm(EventLine2020.shp);rm(EventLine2021.shp);rm(EventLine2022.shp);rm(EventLine2023.shp)
rm(i);rm(path);rm(temp)
gc()

SR <- vect("./Boundary/SouthernRockyBoundary_10kmBuff.shp") ## loading in the study area boundary
SR_FLs <- crop(stacked_FL,SR) ## cropping US fire lines to the study region
gc();rm(stacked_FL)
SR_FLs <- terra::unique(SR_FLs) ## removing duplicate geometries and observations
gc()

vec <- values(SR_FLs) ## making a dataframe to gutcheck the process thus far
colnames(vec) ## Did a double check with all years, confirmed that Feature Cat is only needed column to query
table(vec$DeleteThis)
table(unique(vec$year))
table(vec$FeatureCat)
rm(vec);gc()

SR_FLs <- SR_FLs[SR_FLs$DeleteThis == "No" &
                   SR_FLs$FeatureCat == "Completed Burnout" |
                   SR_FLs$FeatureCat == "Completed Dozer Line" |
                   SR_FLs$FeatureCat == "Completed Fuel Break" | 
                   SR_FLs$FeatureCat == "Completed Hand Line" |
                   SR_FLs$FeatureCat == "Completed Line" |
                   SR_FLs$FeatureCat == "Completed Mixed Construction Line" |
                   SR_FLs$FeatureCat == "Completed Plow Line" |
                   SR_FLs$FeatureCat == "Completed Road as Line" |
                   SR_FLs$FeatureCat == "Contained Line",]
gc()
vec <- values(SR_FLs)
table(vec$FeatureCat)
rm(vec)
gc()

writeVector(SR_FLs, "./NIFC Lines/SR_FLs.shp", overwrite = TRUE) ## writing the shapefile so I can call it later
rm(SR_FLs)

#### Processing Fire Polygon Data from NIFC (2019-2023) ####
temp <- list.files(path = "./NIFC Polygons/",pattern="*.shp") ## creating a vector that has all the files in the working directory with .shp extensions
temp <- temp[1:5] ## setting 1:5 as script will write a compiled .shp file (stored in archived data)

for(i in 1:length(temp)) {
  path <- paste("./NIFC Polygons/", temp[i], sep = "") ## specifying the relative pathway for assign
  assign(temp[i], terra::vect(path)) ## assigning the shapefiles
} ## loading in the shapefiles I want

## cleaning the fire data year by year
## this is done because combining the data prior to cleaning causes issues with computer memory
EventPolygon2019.shp$year <- 2019
length(unique(EventPolygon2019.shp$IncidentNa)) ## 23722 events (whole USA)
EventPolygon2019.shp$IncidentNa <- tolower(gsub("[[:punct:][:space:]]", "", EventPolygon2019.shp$IncidentNa))
EventPolygon2019.shp <- EventPolygon2019.shp[!grepl("rx",EventPolygon2019.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2019.shp <- EventPolygon2019.shp[!grepl("pileburn",EventPolygon2019.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2019.shp <- EventPolygon2019.shp[!grepl("falsealarm",EventPolygon2019.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2019.shp <- EventPolygon2019.shp[!grepl("baer",EventPolygon2019.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2019.shp <- EventPolygon2019.shp[!grepl("test",EventPolygon2019.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2019.shp <- EventPolygon2019.shp[!grepl("delete",EventPolygon2019.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2019.shp <- EventPolygon2019.shp[order(EventPolygon2019.shp$GISAcres, decreasing = TRUE),] ## sorting by size
EventPolygon2019.shp <- EventPolygon2019.shp[!duplicated(EventPolygon2019.shp$IncidentNa),] ## keeping the largest fire per each name (there are errors in here)
length(unique(EventPolygon2019.shp$IncidentNa)) ## 20717
EventPolygon2019.shp <- EventPolygon2019.shp[!EventPolygon2019.shp$GISAcres < 1,]
EventPolygon2019.shp <- EventPolygon2019.shp[EventPolygon2019.shp$DeleteThis == "No",]
length(unique(EventPolygon2019.shp$IncidentNa)) ## 7342
gc()

EventPolygon2020.shp$year <- 2020
length(unique(EventPolygon2020.shp$IncidentNa)) ## 40386
EventPolygon2020.shp$IncidentNa <- tolower(gsub("[[:punct:][:space:]]", "", EventPolygon2020.shp$IncidentNa))
EventPolygon2020.shp <- EventPolygon2020.shp[!grepl("rx",EventPolygon2020.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2020.shp <- EventPolygon2020.shp[!grepl("pileburn",EventPolygon2020.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2020.shp <- EventPolygon2020.shp[!grepl("falsealarm",EventPolygon2020.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2020.shp <- EventPolygon2020.shp[!grepl("baer",EventPolygon2020.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2020.shp <- EventPolygon2020.shp[!grepl("test",EventPolygon2020.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2020.shp <- EventPolygon2020.shp[!grepl("delete",EventPolygon2020.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2020.shp <- EventPolygon2020.shp[order(EventPolygon2020.shp$GISAcres, decreasing = TRUE),] ## sorting by size
EventPolygon2020.shp <- EventPolygon2020.shp[!duplicated(EventPolygon2020.shp$IncidentNa),] ## keeping the largest fire per each name (there are errors in here)
length(unique(EventPolygon2020.shp$IncidentNa)) ## 35882
EventPolygon2020.shp <- EventPolygon2020.shp[!EventPolygon2020.shp$GISAcres < 1,]
EventPolygon2020.shp <- EventPolygon2020.shp[EventPolygon2020.shp$DeleteThis == "No",]
length(unique(EventPolygon2020.shp$IncidentNa)) ## 9517
gc()

EventPolygon2021.shp$year <- 2021
length(unique(EventPolygon2021.shp$IncidentNa)) ## 41647
EventPolygon2021.shp <- EventPolygon2021.shp[EventPolygon2021.shp$IncidentNa != "InTerNaTiOnAl fAlLs WaTeR tOwEr \xed\xa0\xbd\xed\xb7\xbc",] ## this name is causing alot of trouble, so I removed it
EventPolygon2021.shp <- EventPolygon2021.shp[EventPolygon2021.shp$IncidentNa != "fire \xed\xa0\xbd\xed\xb4\xa5 \xed\xa0\xbd\xed\xb3\x9b \xed\xa0\xbd\xed\xb1\xa8‍\xed\xa0\xbd\xed\xba\x92",] ## this name is causing alot of trouble, so I removed it
EventPolygon2021.shp$IncidentNa <- tolower(gsub("[[:punct:][:space:]]", "", EventPolygon2021.shp$IncidentNa))
EventPolygon2021.shp <- EventPolygon2021.shp[!grepl("rx",EventPolygon2021.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2021.shp <- EventPolygon2021.shp[!grepl("pileburn",EventPolygon2021.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2021.shp <- EventPolygon2021.shp[!grepl("falsealarm",EventPolygon2021.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2021.shp <- EventPolygon2021.shp[!grepl("baer",EventPolygon2021.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2021.shp <- EventPolygon2021.shp[!grepl("test",EventPolygon2021.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2021.shp <- EventPolygon2021.shp[!grepl("delete",EventPolygon2021.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2021.shp <- EventPolygon2021.shp[order(EventPolygon2021.shp$GISAcres, decreasing = TRUE),] ## sorting by size
EventPolygon2021.shp <- EventPolygon2021.shp[!duplicated(EventPolygon2021.shp$IncidentNa),] ## keeping the largest fire per each name (there are errors in here)
length(unique(EventPolygon2021.shp$IncidentNa)) ## 37317
EventPolygon2021.shp <- EventPolygon2021.shp[!EventPolygon2021.shp$GISAcres < 1,]
EventPolygon2021.shp <- EventPolygon2021.shp[EventPolygon2021.shp$DeleteThis == "No",]
length(unique(EventPolygon2021.shp$IncidentNa)) ## 15385
gc()

EventPolygon2022.shp$year <- 2022
length(unique(EventPolygon2022.shp$IncidentNa)) ## 33027
EventPolygon2022.shp$IncidentNa <- tolower(gsub("[[:punct:][:space:]]", "", EventPolygon2022.shp$IncidentNa))
EventPolygon2022.shp <- EventPolygon2022.shp[!grepl("rx",EventPolygon2022.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2022.shp <- EventPolygon2022.shp[!grepl("pileburn",EventPolygon2022.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2022.shp <- EventPolygon2022.shp[!grepl("falsealarm",EventPolygon2022.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2022.shp <- EventPolygon2022.shp[!grepl("baer",EventPolygon2022.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2022.shp <- EventPolygon2022.shp[!grepl("test",EventPolygon2022.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2022.shp <- EventPolygon2022.shp[!grepl("delete",EventPolygon2022.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2022.shp <- EventPolygon2022.shp[order(EventPolygon2022.shp$GISAcres, decreasing = TRUE),] ## sorting by size
EventPolygon2022.shp <- EventPolygon2022.shp[!duplicated(EventPolygon2022.shp$IncidentNa),] ## keeping the largest fire per each name (there are errors in here)
length(unique(EventPolygon2022.shp$IncidentNa)) ## 29943
EventPolygon2022.shp <- EventPolygon2022.shp[!EventPolygon2022.shp$GISAcres < 1,]
EventPolygon2022.shp <- EventPolygon2022.shp[EventPolygon2022.shp$DeleteThis == "No",]
length(unique(EventPolygon2022.shp$IncidentNa)) ## 10068
gc()

EventPolygon2023.shp$year <- 2023
length(unique(EventPolygon2023.shp$IncidentNa)) ## 40903
EventPolygon2023.shp$IncidentNa <- tolower(gsub("[[:punct:][:space:]]", "", EventPolygon2023.shp$IncidentNa))
EventPolygon2023.shp <- EventPolygon2023.shp[!grepl("rx",EventPolygon2023.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2023.shp <- EventPolygon2023.shp[!grepl("pileburn",EventPolygon2023.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2023.shp <- EventPolygon2023.shp[!grepl("falsealarm",EventPolygon2023.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2023.shp <- EventPolygon2023.shp[!grepl("baer",EventPolygon2023.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2023.shp <- EventPolygon2023.shp[!grepl("test",EventPolygon2023.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2023.shp <- EventPolygon2023.shp[!grepl("delete",EventPolygon2023.shp$IncidentNa, ignore.case = TRUE),]
EventPolygon2023.shp <- EventPolygon2023.shp[order(EventPolygon2023.shp$GISAcres, decreasing = TRUE),] ## sorting by size
EventPolygon2023.shp <- EventPolygon2023.shp[!duplicated(EventPolygon2023.shp$IncidentNa),] ## keeping the largest fire per each name (there are errors in here)
length(unique(EventPolygon2023.shp$IncidentNa)) ## 33869
EventPolygon2023.shp <- EventPolygon2023.shp[!EventPolygon2023.shp$GISAcres < 1,]
EventPolygon2023.shp <- EventPolygon2023.shp[EventPolygon2023.shp$DeleteThis == "No",]
length(unique(EventPolygon2023.shp$IncidentNa)) ## 10760
gc()

### shortening the number of columns for each of the fire years
EventPolygon2019.shp <- EventPolygon2019.shp[,c(colnames(values(EventPolygon2019.shp)) == "IncidentNa" | 
                                                  colnames(values(EventPolygon2019.shp)) == "GISAcres" |
                                                  colnames(values(EventPolygon2019.shp)) == "year")]
EventPolygon2020.shp <- EventPolygon2020.shp[,c(colnames(values(EventPolygon2020.shp)) == "IncidentNa" | 
                                                  colnames(values(EventPolygon2020.shp)) == "GISAcres" |
                                                  colnames(values(EventPolygon2020.shp)) == "year")]
EventPolygon2021.shp <- EventPolygon2021.shp[,c(colnames(values(EventPolygon2021.shp)) == "IncidentNa" | 
                                                  colnames(values(EventPolygon2021.shp)) == "GISAcres" |
                                                  colnames(values(EventPolygon2021.shp)) == "year")]
EventPolygon2022.shp <- EventPolygon2022.shp[,c(colnames(values(EventPolygon2022.shp)) == "IncidentNa" | 
                                                  colnames(values(EventPolygon2022.shp)) == "GISAcres" |
                                                  colnames(values(EventPolygon2022.shp)) == "year")]
EventPolygon2023.shp <- EventPolygon2023.shp[,c(colnames(values(EventPolygon2023.shp)) == "IncidentNa" | 
                                                  colnames(values(EventPolygon2023.shp)) == "GISAcres" |
                                                  colnames(values(EventPolygon2023.shp)) == "year")]
gc()

## double checking each of the column names
colnames(values(EventPolygon2019.shp))
colnames(values(EventPolygon2020.shp))
colnames(values(EventPolygon2021.shp))
colnames(values(EventPolygon2022.shp))
colnames(values(EventPolygon2023.shp))

head(values(EventPolygon2021.shp)) ## checking the start of the 2021 year to see how it looks
gc()

FirePoly <- rbind(EventPolygon2019.shp,EventPolygon2020.shp,EventPolygon2021.shp,EventPolygon2022.shp,EventPolygon2023.shp)
rm(EventPolygon2019.shp);rm(EventPolygon2020.shp);rm(EventPolygon2021.shp);rm(EventPolygon2022.shp);rm(EventPolygon2023.shp)
gc() ## FirePoly now has all years in the study

table(terra::is.valid(FirePoly)) ## checking the validity of geometery
## have several hundred invalid topologies
FirePoly <- terra::makeValid(FirePoly)
table(terra::is.valid(FirePoly)) ## checking the validity of geometery
## validated
gc()

SR_Fires <- crop(FirePoly,SR)
gc();rm(FirePoly)
SR_Fires <- terra::unique(SR_Fires)
gc()

vec <- values(SR_Fires) ## down to 1285 fires in the Southern Rockies from 2019-2023
colnames(vec)
table(unique(vec$year)) ## double checking we have all the years
table(is.na(vec$GISAcres))  ## double checking there isn't anything funky happening here

plot(SR_Fires) ## checking out the plot

count_vertices <- function(v) {
  n <- nrow(v)
  vertex_counts <- numeric(n)
  for(i in 1:n) {
    coords <- crds(v[i])
    # Subtract 1 because last point repeats first point
    vertex_counts[i] <- nrow(coords) - 1
    progress <- i/n*100
    if (progress %% 10 == 0) {
      print(paste(progress, "% done", sep = ""))
    }
  }
  
  return(vertex_counts)
} ## making a function to count vertices to remove any Rx fires (saved a triangles)

# Remove triangles (3 vertices)
n_vertices <- count_vertices(SR_Fires)
SR_Fires <- SR_Fires[n_vertices > 3,]
length(unique(SR_Fires$IncidentNa))
table(SR_Fires$IncidentNa)
plot(SR_Fires)
gc()

writeVector(SR_Fires,"./NIFC Polygons/SR_Fires.shp", overwrite = TRUE)
rm(SR_Fires);rm(vec);rm(i);rm(n_vertices);rm(path);rm(temp);rm(count_vertices)
gc()

#### Loading Processed Spatial Data for Extraction ####
SR <- vect("./Boundary/SouthernRockyBoundary_10kmBuff.shp") ## loading in the study area boundary, may or may not need to re:load based on where in script you are starting
SR_Fires <- vect("./NIFC Polygons/SR_Fires.shp")
SR_FLs <- vect("./NIFC Lines/SR_FLs.shp")

## Landcover
temp <- list.files(path = "./Predictor Rasters/Land Cover/",pattern="*.tif") ## creating a vector that has all the files in the working directory with .tif extensions
for(i in 1:length(temp)) {
  path <- paste("./Predictor Rasters/Land Cover/", temp[i], sep = "") ## specifying the relative pathway for assign
  assign(temp[i], terra::rast(path)) ## assigning the rasters
} ## loading in the raters I want
LandCov <- do.call(c, lapply(temp, function(x) get(x)))
names(LandCov) <- temp
rm(list = temp)

## Topographic
temp <- list.files(path = "./Predictor Rasters/Topographic/",pattern="*.tif") ## creating a vector that has all the files in the working directory with .tif extensions
for(i in 1:length(temp)) {
  path <- paste("./Predictor Rasters/Topographic/", temp[i], sep = "") ## specifying the relative pathway for assign
  assign(temp[i], terra::rast(path)) ## assigning the rasters
} ## loading in the raters I want
Topo <- do.call(c, lapply(temp, function(x) get(x)))
names(Topo) <- temp
rm(list = temp)

## Climate
temp <- list.files(path = "./Predictor Rasters/Climate/",pattern="*.tif") ## creating a vector that has all the files in the working directory with .tif extensions
for(i in 1:length(temp)) {
  path <- paste("./Predictor Rasters/Climate/", temp[i], sep = "") ## specifying the relative pathway for assign
  assign(temp[i], terra::rast(path)) ## assigning the rasters
} ## loading in the raters I want
Climate <- do.call(c, lapply(temp, function(x) get(x)))
names(Climate) <- temp
rm(list = temp)

rm(i);rm(path);rm(temp)

#### Extracting Environmental Variables for Fire Lines ####
vec <- seq(2019,2023, by = 1) ## change to reflect data range
# i <- 1 ## specifying i because for loop was too memory intensive
Engaged_Lines <- NA

# Pre-allocate list to store results
Engaged_Lines_list <- vector("list", length(vec) * 2)

# Pre-filter data 
Fires_filtered <- SR_Fires[SR_Fires$year %in% vec, ]
FLs_filtered <- SR_FLs[SR_FLs$year %in% vec, ]

## for loop to extract data for fire lines (by year)
for(i in seq_along(vec)){
year_i <- vec[i]

# Filter for current year
Fires_year <- Fires_filtered[Fires_filtered$year == year_i, ]
FLs_year <- FLs_filtered[FLs_filtered$year == year_i, ]

# Buffer operations
Fires_add60 <- buffer(Fires_year, 60) ## adding a buffer around fire perimeter of 60 meters
Fires_minus60 <- buffer(Fires_year, -60)
Fires_EH <- erase(Fires_add60, Fires_minus60)

## getting EH and EF lines
FLs_EH <- terra::intersect(FLs_year, Fires_EH)
FLs_EF <- terra::intersect(FLs_year, Fires_minus60) ## Fires_minus60 represent areas burned over

mid_points_EH <- centroids(FLs_EH) ## finding the midpoints for each line
geometetry_EH <- geom(mid_points_EH) ## extracting the geometery
mid_points_EH <- as.data.frame(mid_points_EH) ## turning into DF for later analysis
perim_EH <- perim(FLs_EH) ## getting the line lengths (in meters)
mid_points_EF <- centroids(FLs_EF) ## same workflow for EF
geometetry_EF <- geom(mid_points_EF)
mid_points_EF <- as.data.frame(mid_points_EF)
perim_EF <- perim(FLs_EF)

# make lines into small polygons
FLs_EH <- buffer(FLs_EH, 1) ## buffering by 1 m to make them small polygons (necessary for exactextractr)
FLs_EF <- buffer(FLs_EF, 1)

# Extract extract with average value for fire lines
# exactextractr returns coverage fraction 
FLs_EH <- tidyterra::as_sf(FLs_EH) ## need to make into sf objects first
FLs_EH <- sf::as_Spatial(FLs_EH$geometry)
FLs_EH_topo <- exact_extract(Topo, FLs_EH, fun = "mean", progress = TRUE)
FLs_EH_cov <- exact_extract(LandCov, FLs_EH, fun = "mean", progress = TRUE)
FLs_EH_clim <- exact_extract(Climate, FLs_EH, fun = "mean", progress = TRUE)

FLs_EF <- tidyterra::as_sf(FLs_EF) ## need to make into sf objects first
FLs_EF <- sf::as_Spatial(FLs_EF$geometry)
FLs_EF_topo <- exact_extract(Topo, FLs_EF, fun = "mean", progress = TRUE)
FLs_EF_cov <- exact_extract(LandCov, FLs_EF, fun = "mean", progress = TRUE)
FLs_EF_clim <- exact_extract(Climate, FLs_EF, fun = "mean", progress = TRUE)

# Combine into dataframes
extracted_FLs_EH <- cbind(mid_points_EH,geometetry_EH,perim_EH,FLs_EH_topo,FLs_EH_cov,FLs_EH_clim)
extracted_FLs_EF <- cbind(mid_points_EF,geometetry_EF,perim_EF,FLs_EF_topo,FLs_EF_cov,FLs_EF_clim)

## getting only the unique geometeries, there are obvious rounding errors if it is not done like this
extracted_FLs_EH$part <- as.integer(base::round(extracted_FLs_EH$x+extracted_FLs_EH$y, digits = 0))
extracted_FLs_EF$part <- as.integer(base::round(extracted_FLs_EF$x+extracted_FLs_EF$y, digits = 0))

EH <- extracted_FLs_EH %>%
  distinct(part, .keep_all = TRUE) ## getting unique observations
EF <- extracted_FLs_EF %>%
  distinct(part, .keep_all = TRUE)

EF$Stat <- "EF"
EH$Stat <- "EH"

colnames(EF) <- colnames(EH) ## syncing colnames for later

Engaged_Lines_list[[2*i - 1]] <- EF
Engaged_Lines_list[[2*i]] <- EH

print(year_i)
gc()
}

## cleaning the extracted data a bit
Engaged_Lines <- do.call(rbind, Engaged_Lines_list);gc()
table(Engaged_Lines$year)
# 2019 2020 2021 2022 2023 
#  236 1994  228 2998  517 
table(Engaged_Lines$Stat)
# EF   EH 
# 1886 4087 
length(unique(Engaged_Lines$IncidentNa.1))
# 56
table(Engaged_Lines$IncidentNa.1) ## a few errors in the fire names

## correcting a few names
Engaged_Lines$IncidentNa.1[Engaged_Lines$IncidentNa.1 == "calfcanyon"] <- "hermitspeak"
Engaged_Lines$IncidentNa.1[Engaged_Lines$IncidentNa.1 == "middleforkfire"] <- "middlefork"
Engaged_Lines$IncidentNa.1[Engaged_Lines$IncidentNa.1 == "nidnight"] <- "midnight"
Engaged_Lines$IncidentNa.1[Engaged_Lines$IncidentNa.1 == "mullenfire"] <- "mullen"
length(unique(Engaged_Lines$IncidentNa.1)) ## 52 total fires

perim_check <- data.frame(name = SR_Fires$IncidentNa,
                          perimeter = perim(SR_Fires)) ## making a df to check line lengths (perim_m)
perim_check$name[perim_check$name == "calfcanyon"] <- "hermitspeak"
perim_check$name[perim_check$name == "middleforkfire"] <- "middlefork"
perim_check$name[perim_check$name == "nidnight"] <- "midnight"
perim_check$name[perim_check$name == "mullenfire"] <- "mullen"

perim_check <- perim_check[order(perim_check$name),]
perim_check <- perim_check[!duplicated(perim_check$name),] ## removing duplicated (lower area observations of same fire)

table(Engaged_Lines$perim_EH >= perim_check$perimeter[match(Engaged_Lines$IncidentNa.1, perim_check$name)])
## this table is telling me there are 4 observations where the line perimeter >= fire perimeter

## looking at other observations
plot(Engaged_Lines$perim_EH)
abline(h = mean(Engaged_Lines$perim_EH) + 3.291*(sd(Engaged_Lines$perim_EH)/(sqrt(length(Engaged_Lines$perim_E))))) ## adding line above the 99.9% confidence intervals (3.291 == t value for CI)
## there are many that are above the 99.9% confidence interval
table(Engaged_Lines$perim_EH > mean(Engaged_Lines$perim_EH) + 3.291*(sd(Engaged_Lines$perim_EH)/(sqrt(length(Engaged_Lines$perim_E))))) ## table of observations above, clearly outliers
## 586 above CI
Engaged_Lines$DeleteThis <- Engaged_Lines$perim_EH > mean(Engaged_Lines$perim_EH) + 3.291*(sd(Engaged_Lines$perim_EH)/(sqrt(length(Engaged_Lines$perim_E)))) ## moving extreme values to "DELETETHIS" columns
table(Engaged_Lines$DeleteThis) ## checking that this matches other table
Engaged_Lines <- Engaged_Lines[Engaged_Lines$DeleteThis == FALSE,] ## only keeping FALSE values

plot(Engaged_Lines$perim_EH) ## this looks more reasonable

## creating a SpatVector that matches the observations of Engaged_Lines
SR_midpoints <- centroids(SR_FLs)
SR_coords <- crds(SR_midpoints)

## Create matching keys
SR_key <- paste(round(SR_coords[,1], 5), round(SR_coords[,2], 5), sep = "_")
Engaged_key <- paste(round(Engaged_Lines$x, 5), 
                     round(Engaged_Lines$y, 5), sep = "_")
stat_lookup <- setNames(Engaged_Lines$Stat, Engaged_key) ## creating a status vector to add to FLs SpatVector

## Subset using %in% 
SR_FLs_trimmed <- SR_FLs[SR_key %in% Engaged_key, ]

## Add the column using the matching keys
SR_FLs_trimmed$stat <- stat_lookup[SR_key[SR_key %in% Engaged_key]]

## saving cleaned fire lines for mapping in ArcPRO
writeVector(SR_FLs_trimmed, "./Results/LinesForMapping.shp")

## formating the create date for matching with VPD, wind, and growth
str(Engaged_Lines$CreateDate)
Engaged_Lines$CreateDate <- as.Date(Engaged_Lines$CreateDate)
plot(Engaged_Lines$CreateDate)
max(Engaged_Lines$CreateDate, na.rm = TRUE)
table(is.na(Engaged_Lines$CreateDate))
Engaged_Lines <- Engaged_Lines[complete.cases(Engaged_Lines$CreateDate),] ## removing the data we have no date for

## which columns to keep 
colnames(Engaged_Lines)
## Stat - 102
## IncidentNa.1 - 39
## year - 31
## XY - 44,45
## CreateDate - 13
## FeatureCat - 2
## perim - 47
## extracted predictors - 48:101

Engaged_Lines <- Engaged_Lines[,c(102,39,31,44,45,13,2,47:101)] ## trimming out the unnecessary columns
colnames(Engaged_Lines) <- c("stat", "fire", "year", "x", "y", "date", "lineconstr", "perim_m","elev",
                             "slope","tpi","aspen_prop","dougfir_prop", "gambel_prop","grass_prop",
                             "lodgepole_prop", "other_prop","pj_prop","pondo_prop", "sf_prop","shrub_prop", 
                             "202303","201904","202004","212104","202204","202304",
                             "201905","202005","022105","202205","202305", 
                             "201906","202006","202106","202206","202306",
                             "201907","202007","202107","202207","202307",
                             "201908","202008","202108","202208","202308",
                             "201909","202009","202109","202209","202309",
                             "201910","202010","202110","202210","202310",
                             "201911","202011","202111","202211","202311") ## renaming VPD columns to match with a date column
colnames(Engaged_Lines)
Engaged_Lines$mdate <- as.character(Engaged_Lines$date) ## making the match date col
Engaged_Lines$mdate <- gsub("-", "",Engaged_Lines$mdate) ## removing dashed lines
Engaged_Lines$mdate <- as.integer(substr(Engaged_Lines$mdate,1,6)) ## pulling out day for matching

Engaged_Lines$VPd <- NA
for(i in 1:nrow(Engaged_Lines)){
  skip_to_next <- FALSE
  tryCatch(Engaged_Lines$VPd[i] <- Engaged_Lines[i,which(colnames(Engaged_Lines) == Engaged_Lines$mdate[i])], error = function(e){skip_to_next <- TRUE})
  if(skip_to_next){next}
  progress <- i/nrow(Engaged_Lines)*100
  if (progress %% 5 == 0) {
    print(paste(progress, "% done", sep = ""))
  }
} ## matching the VPd of the month that line burned 

hist(Engaged_Lines$VPd) ## well that looks right
str(Engaged_Lines)
Engaged_Lines <- as.data.frame(Engaged_Lines)
Engaged_Lines$x <- as.numeric(Engaged_Lines$x)
Engaged_Lines$y <- as.numeric(Engaged_Lines$y)
EL_forTransform <- vect(Engaged_Lines, geom = c("x","y"), crs = crs(SR))
EL_forTransform <- project(EL_forTransform, crs("EPSG:4326"))
Engaged_Lines$lat <- geom(EL_forTransform)[,4]
Engaged_Lines$long <- geom(EL_forTransform)[,3]

## adding wind data
wind <- read.csv("./Predictor Rasters/Wind csv/WindData.csv") ## Reading in Wind Data from Airports in the S. Rockies
wind <- wind[,c(1,3,4,6,7,12)] ## Sub-setting the data to the relevant information: StationID, latitude, longitude, average wind speed for that date (mph), and 5 second maximum wind speed (mph)
colnames(wind) <- c("ID","lat","long","date","avgWind_mph","gustWind_mph") ## renaming variables
wind$avgWind_mph[is.na(wind$avgWind_mph)] <- -9999 ## setting NA values to -9999, as there is a mix of NA and -9999 in the dataframe. Values are reported as -9999 if the datalogger did not record a value during its period of observation. NA values represent missing data outside the period of observation for a datalogger.
wind$gustWind_mph[is.na(wind$gustWind_mph)] <- -9999 ## repeating the same for wind gusts

wind_SNOTEL <- read.csv("./Predictor Rasters/Wind csv/WindData_SNOTEL.csv") ## reading in wind data from the SNOTEL network
wind_SNOTEL <- wind_SNOTEL[,c(2,4,5,1,6,7)]  ## Sub-setting the data to match the 'wind' object
colnames(wind_SNOTEL) <- c("ID","lat","long","date","avgWind_mph","gustWind_mph") ## renaming the columns to match the 'wind' object for a rbind
wind_SNOTEL$avgWind_mph[is.na(wind_SNOTEL$avgWind_mph)] <- -9999 ## setting NA values to -9999 for same reasoning as prior
wind_SNOTEL$gustWind_mph[is.na(wind_SNOTEL$gustWind_mph)] <- -9999

wind_CSU <- read.csv("./Predictor Rasters/Wind csv/WindData_CSUmtnCampus.csv") ## reading in wind data from CSU mountain campus, it needs a little extra processing
wind_CSU <- wind_CSU[,c(4,2,3,1,10,13)] ## subsetting down to the same variables
str(wind_CSU) ## the TIMESTAMP column indicates that wind is reported every 5 minutes, we will need to average by day and get the max values for gusts
wind_CSU$month <- as.numeric(stringr::str_split_fixed(wind_CSU$TIMESTAMP,"/",3)[,1]) ## creating a separate column for the months
wind_CSU$day <- as.numeric(stringr::str_split_fixed(wind_CSU$TIMESTAMP,"/",3)[,2]) ## creating a separate column for the days
wind_CSU$year <- as.numeric(substr(stringr::str_split_fixed(wind_CSU$TIMESTAMP,"/",3)[,3],0,4)) ## creating a separate column for the years

latlong <- c(unique(wind_CSU$lat), unique(wind_CSU$long)) ## making an object that has the lat/long to pair with the data after it is averaged per day
wind_CSU <- wind_CSU %>%
  group_by(year,month,day) %>%
  summarise(avgWind_mph = mean(WindSpeed),
            gustWind_mph = max(HWG_speed)) ## summarizing the the wind_CSU data to create average wind speeds per day and max gusts per day
wind_CSU <- as.data.frame(wind_CSU) ## turning the tibble into a dateframe
wind_CSU$avgWind_mph <- wind_CSU$avgWind_mph*2.237 ## the data was collected in m/s but the rest of the data are in mph, I am converting the units here
wind_CSU$gustWind_mph <- wind_CSU$gustWind_mph*2.237 ## another unit conversion from m/s to mph

wind_RAWS <- read.csv("./Predictor Rasters/Wind csv/RAWS_data_SRockies.csv") ## reading in the RAWS (Remote Automatic Weather Station) data
wind_RAWS <- wind_RAWS[,c(16,17,18,1,6,8)] ## subsetting the data to only be average and gust speeds. We might need furhter data to calculate the energy release component (ERC) but currently do not have the formula
wind_RAWS$Spped.Ave..m.s <- wind_RAWS$Spped.Ave..m.s*2.237 ## data was collected in m/s so I am converting to mph
wind_RAWS$Speed.Gust.m.s <- wind_RAWS$Speed.Gust.m.s*2.237

colnames(wind) ## looking at column names
colnames(wind_CSU) ## seeing the differences so I can rbind in moment
wind_CSU$ID <- "CSU_mtn" ## Adding an ID column
wind_CSU$lat <- latlong[1] ## Adding the lat
wind_CSU$long <- latlong[2] ## Adding the long
wind_CSU$date <- paste(wind_CSU$month,wind_CSU$day,wind_CSU$year, sep = "/") ## repasting and formatting the date information so it is standardized
table(is.na(wind_CSU$gustWind_mph)) ## Double checking that there are no NA values
wind_CSU <- wind_CSU[,c(6,7,8,9,4,5)] ## making sure the columns are in the correct order for the rbind

colnames(wind_RAWS) ## checking the column names
colnames(wind_RAWS) <- c("ID", "lat", "long", "date", "avgWind_mph", "gustWind_mph") ## changing the column names so we can bind
## I might end up using only th RAWS data depending on gaps

wind <- rbind(wind,wind_CSU,wind_SNOTEL,wind_RAWS) ## combining all of the wind data into a singular object

hist(wind$avgWind_mph[wind$avgWind_mph >= 0]) ## looking at a distribution of the data (excluding weird negative values)
length(wind$avgWind_mph[wind$avgWind_mph < 0]) ## non-usable values
hist(wind$gustWind_mph[wind$gustWind_mph >= 0])
length(wind$gustWind_mph[wind$gustWind_mph < 0]) ## also non-usable values
table(is.na(wind$avgWind_mph)) ## no NA values for each
table(is.na(wind$gustWind_mph)) ## looking at which stations are responsible for the NA values (mostly BARTLEY in NM)

wind$avgWind_mph[wind$avgWind_mph < 0] <- -9999 ## turning negative values to -9999 then later NA as there is no negative wind speed
wind$gustWind_mph[wind$gustWind_mph <0] <- -9999

length(unique(wind$ID)) ## We have wind data for 101 different stations from RAWS
## plotting the distribution of wind stations and sampled firelines
plot(unique(wind$long), unique(wind$lat), col = "red", pch = 16, cex = 0.5,
     xlab = "longitude",
     ylab = "latitude")
points(Engaged_Lines$long,Engaged_Lines$lat, col = "black", pch = 16, cex = 0.5)
points(unique(wind$long), unique(wind$lat), col = "red", pch = 16, cex = 0.5)
legend("topleft",legend = c("stations", "sampled points"), pch = 16, col = c("red","black"), bty = "n")

## nearest neighbor analysis
## I am matching the sampled fire line points to the nearest wind station and attributing the wind value from the station on the day that the fire line was engaged to the sampled point
stations <- matrix(c(wind$lat[match(unique(wind$ID),wind$ID)],wind$long[match(unique(wind$ID),wind$ID)]),nrow = length(unique(wind$ID)),ncol = 2, byrow = FALSE ) ## creating a matrix of lat/long for the wind stations
rownames(stations) <- unique(wind$ID) ## naming the rows for the wind station it is associated with, this is useful for matching later on
firepoints  <- matrix(c(Engaged_Lines$lat,Engaged_Lines$long),nrow = nrow(Engaged_Lines),ncol = 2, byrow = FALSE) ## making a matrix of the points for each fire line so I look at the distance between a sample point and each station
m1 <- stations[,c(2,1)] ## creating a second matrix with a shorter name that is long/lat (as opposed to lat/long) for the wind stations
m2 <- firepoints[,c(2,1)] ## creating a second matrix with a shorter name that is long/lat (as opposed to lat/long) for the sampled fire lines
dist.mat <- distance(m1,m2, lonlat = TRUE) ## the distance function (terra) calculates distance between points
rownames(dist.mat) <- rownames(stations) ## naming the rows so I can determine the closest station to each point
dist.df <- as.data.frame(dist.mat) ## turning the matrix into a data frame so I can use apply functions
dists <- as.vector(apply(dist.df,2,min)) ## creating a vector that has the shortest distance to a station for each sampled fire line
max(dists)/1000 ## furthest away a station is from a point 83 km (using the midpoint of a line)
min(dists)/1000 ## closest is 0.4 km (using the midpoint of a line)

hist(dists/1000,
     main = "Histogram of Distance to Weather Station",
     las = 1,
     xlab = "Distance (km)") ## histogram of the distances

closest.stations <- NA ## creating an empty vector that I can assign a closest station for each poin
for(i in 1:ncol(dist.mat)){
  closest.stations[i] <- rownames(dist.mat)[!is.na(match(dist.mat[,i],dists))]
} ## identifying the station with the smallest distance to each sampled fire lines point

length(unique(closest.stations))  ## 41 stations were used (out of the 101 available)
table(is.na(closest.stations)) ## no NA values
Engaged_Lines$Station <- closest.stations ## Assigning each point its closest station that has wind data

## Configuring the dates column in the wind object to facilitate searching
table(is.na(wind$date)) ## no NA values
wind$date <- paste(as.numeric(stringr::str_split_fixed(wind$date,"/",3)[,3]),
                   ifelse(as.numeric(stringr::str_split_fixed(wind$date,"/",3)[,1]) < 10 ,paste(0,as.numeric(stringr::str_split_fixed(wind$date,"/",3)[,1]),sep=""),as.numeric(stringr::str_split_fixed(wind$date,"/",3)[,1])),
                   ifelse(as.numeric(stringr::str_split_fixed(wind$date,"/",3)[,2]) < 10 ,paste(0,as.numeric(stringr::str_split_fixed(wind$date,"/",3)[,2]),sep=""),as.numeric(stringr::str_split_fixed(wind$date,"/",3)[,2])),
                   sep = "") ## updating the date to be YYYYMMDD format instead of MM/DD/YYYY

Engaged_Lines$mdate <- as.character(Engaged_Lines$date) ## making the match date col
Engaged_Lines$mdate <- as.integer(gsub("-", "",Engaged_Lines$mdate)) ## removing dashed lines

## Matching the stations and days to individual points to give each point wind data!
Engaged_Lines$avgWind_mph <- NA ## creating blank columns to add data to from for loop for avg wind speed
Engaged_Lines$gustWind_mph <- NA ## creating blank columns to add data to from for loop for wind gust speed

for(i in 1:nrow(Engaged_Lines)){
  skip_to_next <- FALSE
  tryCatch(Engaged_Lines$avgWind_mph[i] <- wind$avgWind_mph[which(wind$date == Engaged_Lines$mdate[i] & wind$ID == Engaged_Lines$Station[i])], error = function(e) {skipe_to_next <- TRUE})
  tryCatch(Engaged_Lines$gustWind_mph[i] <- wind$gustWind_mph[which(wind$date == Engaged_Lines$mdate[i] & wind$ID == Engaged_Lines$Station[i])], error = function(e) {skipe_to_next <- TRUE})
  if(skip_to_next) {next}
} ## for loop that assigns a data based on matching station ID and date, there is a skip to next function in case there are NA values in the wind data

table(is.na(Engaged_Lines$avgWind_mph)) ## 29 needed to be skipped over
table(is.na(Engaged_Lines$gustWind_mph)) ## 29 skipped for this variable too

table(Engaged_Lines$avgWind_mph == -9999 | is.na(Engaged_Lines$avgWind_mph)) ## 4657 observations with usable data for average wind speeds (using the 1 km point density for fire line sampling points)
table(Engaged_Lines$gustWind_mph == -9999 | is.na(Engaged_Lines$gustWind_mph)) ## 4426 observation with usable data for gust speeds (using the 1 km point density for fire line sampling points)

Engaged_Lines$avgWind_mph[Engaged_Lines$avgWind_mph == -9999] <- NA ## setting -9999 values back to NA for the upcoming modelling
Engaged_Lines$gustWind_mph[Engaged_Lines$gustWind_mph == -9999] <- NA ## setting -9999 values back to NA for the upcoming modelling

hist(Engaged_Lines$avgWind_mph)
hist(Engaged_Lines$gustWind_mph)

hist(Engaged_Lines$gustWind_mph[Engaged_Lines$stat == "EH" & Engaged_Lines$aspen_prop >= 0.5],
     col = rgb(1,1,0,0.5),
     xlim = c(0,100),
     main = "",
     las = 1)
hist(Engaged_Lines$gustWind_mph[Engaged_Lines$stat == "EF" & Engaged_Lines$aspen_prop >= 0.5],
     col = rgb(0,0,1,0.5),
     add = T)

## adding fire growth variable
SIT209 <- read.csv("./Predictor Rasters/ICS-209 csv/ics209-plus-wf_sitreps_1999to2023.csv")
SIT209 <- SIT209[SIT209$CY == 2019 | SIT209$CY == 2020 | SIT209$CY == 2021 | SIT209$CY == 2022 | SIT209$CY == 2023,] ## pulling out the years of interest
colnames(SIT209)
which(colnames(SIT209) == "REPORT_FROM_DATE")
which(colnames(SIT209) == "REPORT_TO_DATE")
which(colnames(SIT209) == "CURR_INCIDENT_AREA")
which(colnames(SIT209) == "INCIDENT_NAME")
SIT209 <- SIT209[,c(106,107,13,54)] ## pulling out the columns of interest
SIT209$fmatch <- paste(gsub(" ","",tolower(SIT209$INCIDENT_NAME)), as.numeric(substr(stringr::str_split_fixed(SIT209$REPORT_FROM_DATE,"/",3)[,1],0,4)), sep = "") 
## making a column to match the fire information between dataframe that includes the names of the fires and the years that they occurred
Engaged_Lines$fmatch <- paste(Engaged_Lines$fire, Engaged_Lines$year, sep = "")

## pulling out the matches from the SIT209 data
f1 <- SIT209[!is.na(match(SIT209$fmatch, Engaged_Lines$fmatch)),] ## pulling out all of the observations from the ICS-209 data that isn't NA for fire events that are present in the FireData dataframe
length(unique(f1$fmatch)) ## there are 40 fires with data in the ICS-209 form that are found in the Engaged_Lines dataframe

f1$month <- as.numeric(stringr::str_split_fixed(f1$REPORT_FROM_DATE,"-",3)[,2])
f1$day <- as.numeric(substr(stringr::str_split_fixed(f1$REPORT_FROM_DATE,"-",3)[,3],0,2))
f1$year <- as.numeric(stringr::str_split_fixed(f1$REPORT_FROM_DATE,"-",3)[,1])  ## I am breaking apart the date column to reorder the dataframe more easily

f2 <- f1 %>% group_by(fmatch, year, month, day) %>%
  summarize(area = sum(CURR_INCIDENT_AREA)) ## this code sums the total incident area per day and organizes the data by event, then year, month, and day
f2 <- as.data.frame(f2) ## turning the tibble back into a dataframe

par(mfrow = c(1,1))
vec <- unique(f2$fmatch)
for(i in 1:length(vec)){
  plot(f2$area[f2$fmatch == vec[i]],
       main = vec[i],
       ylab = "area")
  print(unique(f2$fmatch[f2$fmatch == vec[i]]))
  print(max(f2$area[f2$fmatch == vec[i]]))
} ## printing out the name of each fire and the area
## Go back and look through each figure to assess errors, there are several
## manually resolving errors in the data (there has to be a better way, but I don't know it)
## errors are generally resolved by taking the average between the point preceeding and following it
## if there was ever a question, I googled the fire

## 403 
f2$area[f2$fmatch == "4032023"][4] <- mean(f2$area[f2$fmatch == "4032023"][c(3,5)]) 

## american mesa
f2$area[f2$fmatch == "americanmesa2023"][9] <- mean(f2$area[f2$fmatch == "americanmesa2023"][c(8,10)]) 

## black feather
f2$area[f2$fmatch == "blackfeather2023"][2]  <- mean(f2$area[f2$fmatch == "blackfeather2023"][c(1,3)])

## calwood
f2$area[f2$fmatch == "calwood2020"][3] <- 8960
f2$area[f2$fmatch == "calwood2020"][20] <- f2$area[f2$fmatch == "calwood2020"][19]

## cameron peak
f2$area[f2$fmatch == "cameronpeak2020"][19] <- mean(f2$area[f2$fmatch == "cameronpeak2020"][c(18,20)])
f2$area[f2$fmatch == "cameronpeak2020"][20:40] ## just checking which number the errors are
f2$area[f2$fmatch == "cameronpeak2020"][27] <- f2$area[f2$fmatch == "cameronpeak2020"][26]
f2$area[f2$fmatch == "cameronpeak2020"][34] <- f2$area[f2$fmatch == "cameronpeak2020"][33]
f2$area[f2$fmatch == "cameronpeak2020"][38] <- mean(f2$area[f2$fmatch == "cameronpeak2020"][c(37,39)])
f2$area[f2$fmatch == "cameronpeak2020"][58] <- mean(f2$area[f2$fmatch == "cameronpeak2020"][c(57,59)])

## cerro pelado
f2$area[f2$fmatch == "cerropelado2022"][2]<- mean(f2$area[f2$fmatch == "cerropelado2022"][c(1,3)])
f2$area[f2$fmatch == "cerropelado2022"][7]<- mean(f2$area[f2$fmatch == "cerropelado2022"][c(6,8)])
f2$area[f2$fmatch == "cerropelado2022"][14] <- mean(f2$area[f2$fmatch == "cerropelado2022"][c(13,15)])
f2$area[f2$fmatch == "cerropelado2022"][22] <- mean(f2$area[f2$fmatch == "cerropelado2022"][c(21,23)])
f2$area[f2$fmatch == "cerropelado2022"][51] <- f2$area[f2$fmatch == "cerropelado2022"][50]
f2$area[f2$fmatch == "cerropelado2022"][59] <- f2$area[f2$fmatch == "cerropelado2022"][60]

# ## cooks peak 
f2$area[f2$fmatch == "cookspeak2022"][7] <- mean(f2$area[f2$fmatch == "cookspeak2022"][c(6,8)])

## decker
f2$area[f2$fmatch == "decker2019"][28] <- mean(f2$area[f2$fmatch == "decker2019"][c(27,29)])
f2$area[f2$fmatch == "decker2019"][31] <- mean(f2$area[f2$fmatch == "decker2019"][c(30,32)])
f2$area[f2$fmatch == "decker2019"][49] <- f2$area[f2$fmatch == "decker2019"][50]

## east canyon
## this one is particularly off. I will make it follow 'normal' fire behavior; only burned 2905 acres total
f2$area[f2$fmatch == "eastcanyon2020" & f2$area > 3000] <- 2905

## east troublesome
f2$area[f2$fmatch == "easttroublesome2020"][23] <- f2$area[f2$fmatch == "easttroublesome2020"][24]

## el valle
f2$area[f2$fmatch == "elvalle2023"][1] <- mean(c(f2$area[f2$fmatch == "elvalle2023"][2],1))
f2$area[f2$fmatch == "elvalle2023"][c(3,8)] <- f2$area[f2$fmatch == "elvalle2023"][9]

## grizzly creek
f2$area[f2$fmatch == "grizzlycreek2020"][2] <- mean(f2$area[f2$fmatch == "grizzlycreek2020"][c(1,3)])
f2$area[f2$fmatch == "grizzlycreek2020"][5] <- mean(f2$area[f2$fmatch == "grizzlycreek2020"][c(4,6)])
f2$area[f2$fmatch == "grizzlycreek2020"][7] <- mean(f2$area[f2$fmatch == "grizzlycreek2020"][c(6,8)])
f2$area[f2$fmatch == "grizzlycreek2020"][14] <- mean(f2$area[f2$fmatch == "grizzlycreek2020"][c(13,15)])
f2$area[f2$fmatch == "grizzlycreek2020"][21] <- f2$area[f2$fmatch == "grizzlycreek2020"][22]
f2$area[f2$fmatch == "grizzlycreek2020"][36] <- f2$area[f2$fmatch == "grizzlycreek2020"][37]

## hermits peak
f2$area[f2$fmatch == "hermitspeak2022"][63] <- mean(f2$area[f2$fmatch == "hermitspeak2022"][c(62,64)])
f2$area[f2$fmatch == "hermitspeak2022"][142] <- f2$area[f2$fmatch == "hermitspeak2022"][141]

## high park
f2$area[f2$fmatch == "highpark2022"][2] <- mean(f2$area[f2$fmatch == "highpark2022"][c(1,3)])

## lefthand 
f2$area[f2$fmatch == "lefthand2020"][2] <- mean(f2$area[f2$fmatch == "lefthand2020"][c(1,3)])

## little mesa
f2$area[f2$fmatch == "littlemesa2023"][3] <- mean(f2$area[f2$fmatch == "littlemesa2023"][c(2,4)])
f2$area[f2$fmatch == "littlemesa2023"][5] <- mean(f2$area[f2$fmatch == "littlemesa2023"][c(4,6)])

## lowline
f2$area[f2$fmatch == "lowline2023"][1] <- mean(c(f2$area[f2$fmatch == "lowline2023"][2],1))
f2$area[f2$fmatch == "lowline2023"][11] <- mean(f2$area[f2$fmatch == "lowline2023"][c(10,12)])
f2$area[f2$fmatch == "lowline2023"][c(16:18)] <- f2$area[f2$fmatch == "lowline2023"][19]

# ## luna
f2$area[f2$fmatch == "luna2020"][2] <- mean(f2$area[f2$fmatch == "luna2020"][c(1,3)])
f2$area[f2$fmatch == "luna2020"][4] <- mean(f2$area[f2$fmatch == "luna2020"][c(3,5)])

## middle fork
f2$area[f2$fmatch == "middlefork2020"][34] <- mean(f2$area[f2$fmatch == "middlefork2020"][c(33,35)])
f2$area[f2$fmatch == "middlefork2020"][59] <- f2$area[f2$fmatch == "middlefork2020"][60]

## middle mamm
f2$area[f2$fmatch == "middlemamm2019"][5] <- mean(f2$area[f2$fmatch == "middlemamm2019"][c(4,6)])
f2$area[f2$fmatch == "middlemamm2019"][13] <- mean(f2$area[f2$fmatch == "middlemamm2019"][c(12,15)])
f2$area[f2$fmatch == "middlemamm2019"][14] <- mean(f2$area[f2$fmatch == "middlemamm2019"][c(13,15)])
f2$area[f2$fmatch == "middlemamm2019"][17] <- mean(f2$area[f2$fmatch == "middlemamm2019"][c(16,18)])
f2$area[f2$fmatch == "middlemamm2019"][22] <- mean(f2$area[f2$fmatch == "middlemamm2019"][c(21,23)])

## monday creek
f2$area[f2$fmatch == "mondaycreek2022"][1] <- f2$area[f2$fmatch == "mondaycreek2022"][2]/2

## morgan creek
f2$area[f2$fmatch == "morgancreek2021"][1] <- mean(c(f2$area[f2$fmatch == "morgancreek2021"][2],1))
f2$area[f2$fmatch == "morgancreek2021"][3] <- mean(f2$area[f2$fmatch == "morgancreek2021"][c(2,4)])
f2$area[f2$fmatch == "morgancreek2021"][6] <- mean(f2$area[f2$fmatch == "morgancreek2021"][c(5,7)])
f2$area[f2$fmatch == "morgancreek2021"][11] <- mean(f2$area[f2$fmatch == "morgancreek2021"][c(10,12)])
f2$area[f2$fmatch == "morgancreek2021"][41] <- f2$area[f2$fmatch == "morgancreek2021"][40]

## mullen
f2$area[f2$fmatch == "mullen2020"][8] <- mean(f2$area[f2$fmatch == "mullen2020"][c(7,9)])
f2$area[f2$fmatch == "mullen2020"][44:45] <- f2$area[f2$fmatch == "mullen2020"][46]
f2$area[f2$fmatch == "mullen2020"][51:52] <- f2$area[f2$fmatch == "mullen2020"][53]
f2$area[f2$fmatch == "mullen2020"][55] <- f2$area[f2$fmatch == "mullen2020"][54]

## plumtaw
f2$area[f2$fmatch == "plumtaw2022"][1] <- f2$area[f2$fmatch == "plumtaw2022"][2]/2
f2$area[f2$fmatch == "plumtaw2022"][10:14] <- f2$area[f2$fmatch == "plumtaw2022"][14]

## reveille
f2$area[f2$fmatch == "reveille2019"][3] <- mean(f2$area[f2$fmatch == "reveille2019"][c(2,4)])
f2$area[f2$fmatch == "reveille2019"][16] <- mean(f2$area[f2$fmatch == "reveille2019"][c(15,17)])
f2$area[f2$fmatch == "reveille2019"][22] <- mean(f2$area[f2$fmatch == "reveille2019"][c(21,23)])
f2$area[f2$fmatch == "reveille2019"][24] <- mean(f2$area[f2$fmatch == "reveille2019"][c(23,25)])
f2$area[f2$fmatch == "reveille2019"][33] <- mean(f2$area[f2$fmatch == "reveille2019"][c(32,34)])
f2$area[f2$fmatch == "reveille2019"][c(46,47)] <- f2$area[f2$fmatch == "reveille2019"][48]

## simms
f2$area[f2$fmatch == "simms2022"][1] <- f2$area[f2$fmatch == "simms2022"][2] 
f2$area[f2$fmatch == "simms2022"][5] <- mean(f2$area[f2$fmatch == "simms2022"][c(4,6)])

## spring creek
f2$area[f2$fmatch == "springcreek2023"][5] <- f2$area[f2$fmatch == "springcreek2023"][6]
f2$area[f2$fmatch == "springcreek2023"][19] <- f2$area[f2$fmatch == "springcreek2023"][20]
f2$area[f2$fmatch == "springcreek2023"][c(54,56,57,60,71)] <- f2$area[f2$fmatch == "springcreek2023"][61]

## sylvan
f2$area[f2$fmatch == "sylvan2021"][2] <- mean(f2$area[f2$fmatch == "sylvan2021"][c(1,3)])
f2$area[f2$fmatch == "sylvan2021"][c(20,30)] <- f2$area[f2$fmatch == "sylvan2021"][21]

## trail springs
f2$area[f2$fmatch == "trailsprings2023"][19] <- f2$area[f2$fmatch == "trailsprings2023"][18]

## williams fork
f2$area[f2$fmatch == "williamsfork2020" & f2$month == 8 & f2$day == 15] <- 3963.5

## ymca
f2$area[f2$fmatch == "ymca2020"][4:6] <- mean(f2$area[f2$fmatch == "ymca2020"][c(3,7)])

par(mfrow = c(1,1))
vec <- unique(f2$fmatch)
for(i in 1:length(vec)){
  plot(f2$area[f2$fmatch == vec[i]],
       main = vec[i],
       ylab = "area")
  print(unique(f2$fmatch[f2$fmatch == vec[i]]))
  print(max(f2$area[f2$fmatch == vec[i]]))
} ## printing out the name of each fire and the area
## going back and double checking that they all look okay

vec <- unique(f2$fmatch) ## creating a vector to run the for loop for calculating the daily growth for each unique fire event
for(i in 1:length(vec)){
  tmp <- f2[f2$fmatch == vec[i],]
  tmp$growth <- tmp$area - c(0,tmp$area[c(1:(length(tmp$area))-1)])
  f2$growth[f2$fmatch == vec[i]] <- tmp$growth
} ## for loop to calculate growth. The formula for growth is difference in area at time t from the area at time t+1. For day 1, the area area in time t == 0
## positive values indicate increases in area
## negative values indicate decreases in area
## 0 values indicate no change in area
f2$date <- paste(f2$year,ifelse(f2$month < 10, paste("0", as.character(f2$month), sep = ""), as.character(f2$month)),ifelse(f2$day < 10, paste("0", as.character(f2$day), sep = ""), as.character(f2$day)), sep = "") ## creating a date column that matches the YYYYMMDD format of the mdate column
FireGrowth <- f2[,c(1,7,6)] ## creating a dataframe for only the relevant information pertaining to FireGrowth: FireName, Date, Area growth on that date for that fire
hist(FireGrowth$growth)

## Matching the fire growth for each event on each specific date
vec <- unique(Engaged_Lines$fmatch) ## making a vector for unique event
Engaged_Lines$Growth <- NA ## creating an empty column in the Engaged_Lines dataframe that I can then store the growth covariate from the FireGrowth dataframe
for(i in 1:length(vec)){
  tmp <- Engaged_Lines[Engaged_Lines$fmatch == vec[i],]
  tmp2 <- FireGrowth[FireGrowth$fmatch == vec[i],]
  Engaged_Lines$Growth[Engaged_Lines$fmatch == vec[i]] <- tmp2$growth[(match(tmp$mdate, tmp2$date))]
} ## this code will work with additional fire years with additional ICS-209 forms. For now there are a lot of NA values
hist(Engaged_Lines$Growth) ## there are some extreme events but mostly moderate growth
table(is.na(Engaged_Lines$Growth)) ## 134 observations have no growth data
table(Engaged_Lines$fire[is.na(Engaged_Lines$Growth)]) ## split across a few fires

FireData <- Engaged_Lines[complete.cases(Engaged_Lines),]
table(FireData$stat) ## Sample Size of 4295
## 1482 failed
## 2813 held
table(FireData$stat)/nrow(FireData)*100 ## a little unbalanced
# EF       EH 
# 34.50524 65.49476
length(unique(FireData$fire)) ## 36 total fires

FireData$Jdate <- format(strptime(FireData$mdate, format = "%Y%m%d"), "%j")
FireData$Jdate <- as.integer(FireData$Jdate)

FireData$fire[FireData$Growth == 115512]
FireData$Growth[FireData$Growth == 115512] <- 0 ## removing Growth outlier (real growth but in opposite direction of firelines)

## cleaning up the FireData
colnames(FireData)
## 1-21,64:66,68,69,71,72
FireData <- FireData[,c(1:21,64:66,68,69,71,72)]
write.csv(FireData, "./Results/CleanedFireLines.csv")
rm(list = ls()) # cleans entire global environment
gc()

#### Visualizing Cleaned Fire Lines Data ####
## summary of data
FireData <- read.csv("./Results/CleanedFireLines.csv")

length(unique(FireData$fire))
## 36 fires
table(FireData$stat)
# EF   EH 
# 1482 2813 
(table(FireData$stat)/nrow(FireData))*100
## 35/65 split
length(unique(FireData$fire))
table(is.na(FireData$stat))

groupedCounts <- as.data.frame(FireData %>% group_by(fire) %>% count(stat)) ## getting the count of failed/held lines by fire
EF <- groupedCounts[groupedCounts$stat == "EF",] ## making a separate dataframe for EH/EF lines
EH <- groupedCounts[groupedCounts$stat == "EH",]
length(unique(FireData$fire)[!(unique(FireData$fire) %in% EF$fire)]) ## 9 fires don't have failed lines
length(unique(FireData$fire)[!(unique(FireData$fire) %in% EH$fire)]) ## 1 fires doesn't have held lines
EF[28:36,] <- NA ## adding in extra rows, to bind in the missing fires
EH[36,] <- NA
EF$fire[is.na(EF$fire)] <- unique(FireData$fire)[!(unique(FireData$fire) %in% EF$fire)] ## adding in the missing fires
EH$fire[is.na(EH$fire)] <- unique(FireData$fire)[!(unique(FireData$fire) %in% EH$fire)] 
EF$n[is.na(EF$n)] <- 0 ## adding the count of na files to the EH/EF data.frames
EH$n[is.na(EH$n)] <- 0

summaryTable <- data.frame(fire = unique(FireData$fire),
                           year = FireData$year[duplicated(FireData$fire) == FALSE],
                           EH = EH$n[match(unique(FireData$fire),EH$fire)] , 
                           EF = EF$n[match(unique(FireData$fire),EF$fire)] )
summaryTable <- summaryTable[order(summaryTable$fire),]
sum(summaryTable$EH)
sum(summaryTable$EF)

SR_Fires <- vect("./NIFC Polygons/SR_Fires.shp")
SR_Fires <- as.data.frame(SR_Fires)
summaryTable$area_ha <- SR_Fires$GISAcres[match(summaryTable$fire, SR_Fires$IncidentNa)]
summaryTable$area_ha <- summaryTable$area_ha/2.471 ## acres to ha
summaryTable$area_ha <- round(summaryTable$area_ha,0) ## rounding
hist(summaryTable$area_ha) ## one obvious error
max(summaryTable$area_ha) ## 700k ha is not the correct size
summaryTable$fire[which(summaryTable$area_ha == max(summaryTable$area_ha))] ## cameronpeak is incorrect
## fire size was actually 84544 ha
summaryTable$area_ha[summaryTable$fire == "cameronpeak"] <- 84544
table(summaryTable$area_ha) ## next largest should be hermit peak
summaryTable$fire[which(summaryTable$area_ha == max(summaryTable$area_ha))] 
hist(summaryTable$area_ha)
## rest look approx. correct
summaryTable <- summaryTable[,c(1,5,2,3,4)]

avg_perim <-  stats::aggregate(FireData$perim_m ~ FireData$fire, FUN = mean)[2] ## writing it as a vector first
avg_perim <- as.vector(round(avg_perim$`FireData$perim_m`,0)) ## this avoids a weird column name and odd structure after writing to csv
summaryTable$avg_linelength_m <- avg_perim
write.csv(summaryTable, "./Results/Table1.csv")

sum(summaryTable$EH)/nrow(FireData)*100
## 65% of lines held
vec <- c("cameronpeak","mullen","easttroublesome","hermitspeak") ## 4 largest fires
subset.FD <- FireData[!FireData$fire %in% vec,] ## subsetting out megafires
subset.ST <- summaryTable[!summaryTable$fire %in% vec,]
sum(subset.ST$EH)/nrow(subset.FD)*100
## 82% of lines held when not considering 4 largest fires
round(sum(summaryTable$EF)/sum(subset.ST$EF),0)
## a 6 fold increase in the number of line failures when considering all fires vs. subset

str(summaryTable)
tots <- data.frame(year = c(2019,2020,2021,2022,2023),
                   n.fires = c(length(summaryTable$area_ha[summaryTable$year == 2019]),
                               length(summaryTable$area_ha[summaryTable$year == 2020]),
                               length(summaryTable$area_ha[summaryTable$year == 2021]),
                               length(summaryTable$area_ha[summaryTable$year == 2022]),
                               length(summaryTable$area_ha[summaryTable$year == 2023])),
                   tot.area = c(sum(summaryTable$area_ha[summaryTable$year == 2019]),
                                sum(summaryTable$area_ha[summaryTable$year == 2020]),
                                sum(summaryTable$area_ha[summaryTable$year == 2021]),
                                sum(summaryTable$area_ha[summaryTable$year == 2022]),
                                sum(summaryTable$area_ha[summaryTable$year == 2023])),
                   tot.EF = c(sum(summaryTable$EF[summaryTable$year == 2019]),
                              sum(summaryTable$EF[summaryTable$year == 2020]),
                              sum(summaryTable$EF[summaryTable$year == 2021]),
                              sum(summaryTable$EF[summaryTable$year == 2022]),
                              sum(summaryTable$EF[summaryTable$year == 2023])),
                   tot.EH = c(sum(summaryTable$EH[summaryTable$year == 2019]),
                              sum(summaryTable$EH[summaryTable$year == 2020]),
                              sum(summaryTable$EH[summaryTable$year == 2021]),
                              sum(summaryTable$EH[summaryTable$year == 2022]),
                              sum(summaryTable$EH[summaryTable$year == 2023])))
par(oma = c(0,1.5,0,1.5))
barplot(height = tots$tot.area/1000, names = tots$year,
        las = 1,
        ylim = c(0,300),
        xlab = "")
mtext("Total Area Burned (1,000 ha)", side = 2, line = 4)
text(x = c(0.75,1.9,3.1,4.3,5.5),
     y = (tots$tot.area/1000 + 15),
     cex = 1,
     labels = tots$n.fires)

fl.mat <- as.matrix(tots[,c(4,5)])
rownames(fl.mat) <- c("2019", "2020", "2021","2022", "2023")
colnames(fl.mat) <- c("EF", "EH")

par(oma = c(0,0,0,0))
barplot(t(fl.mat), beside = T,
        las = 2,
        ylim = c(0,1600),
        xlab = "year",
        ylab = "number of lines",
        col = c("darkgoldenrod","navy"))
legend("topright",legend = c("Failed Lines", "Held Lines"), pch = 15, col = c("darkgoldenrod","navy"), bty = "n")
text(x = c(1.5,4.5,7.5,10.4,13.5,2.5,5.5,8.5,11.5,14.5),
     y = (as.vector(fl.mat) + 100),
     cex = 1,
     labels = as.vector(fl.mat))

rm(fl.mat);rm(SR_Fires);rm(summaryTable);rm(tots);rm(i);rm(vec)
rm(EF);rm(EH);rm(groupedCounts);rm(subset.FD);rm(subset.ST);rm(avg_perim)

#### Variable Checking ####
colnames(FireData)
preds <- as.matrix(FireData[,c(9:23,26:28)])
M <- cor(preds, method = c("spearman"))
correlated <- as.numeric(M[which(M > 0.75)])
correlated <- correlated[correlated != 1] 
which(M == correlated, arr.ind = TRUE)
## wind is correlated with wind
M[which(M < -0.75)]

rm(preds);rm(M);rm(correlated)

## Checking for Collinearity
## Custum function from Zuur for collineariy
# Library files for courses provided by: Highland Statistics Ltd.
# To cite these functions, use:
# Mixed effects models and extensions in ecology with R. (2009).
# Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.
# Copyright Highland Statistics LTD.
# VIF FUNCTION.
# To use:  corvif(YourDataFile)
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}


#Support function for corvif. Will not be called by the user
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}
#END VIF FUNCTIONS

colnames(FireData) ## colnames of interst - I am grabbing the number instead
# "fire"           "year"           "lineconstr"     "perim_m"       
# "elev"           "slope"          "tpi"            "aspen_prop"    
# "dougfir_prop"   "gambel_prop"    "grass_prop"     "lodgepole_prop"
# "other_prop"     "pj_prop"        "pondo_prop"     "sf_prop"       
# "shrub_prop"     "VPd"            "avgWind_mph"    "gustWind_mph"  
# "Growth"
MyVar <- colnames(FireData)[c(3,4,8:23,26:28)]
MyVar
corvif(FireData[,MyVar]) ## FireName explains a lot of the variation, this makes makes since given the heirarchical nature of the data
rm(MyVar);rm(corvif);rm(myvif)

#### Frequentest Models ####
pseudo.R.squared <- function(model){
  1 - (model$deviance/model$null.deviance)
} ## R squared function 

## Model 1 - simplest model no random effects
colnames(FireData)
FDmod <- FireData[,c(2,3,13,28,23)]
FDmod$stat <- as.integer(ifelse(FDmod$stat == "EF", 0,1))
m1 <- glm(stat ~ aspen_prop+Growth+VPd,
          family = binomial(link = "logit"),
          data = FDmod)
summary(m1) ## AIC 5832.9

plot(residuals(m1),
     ylab = "residuals",
     pch = 16,
     cex = 0.75,
     col = rgb(0,0,0,0.75),
     las = 1,
     xlab = "index") ## residuals real pretty bad
abline(h = 0, col = "red", lwd = 2, lty = 2)
qqnorm(residuals(m1), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(m1), col = "red", lwd = 2, lty = 2) ## wow horrible, who would have thought /s

## Model 2 - simplest model with random effects
m2 <- lme4::glmer(stat ~scale(aspen_prop)+scale(Growth)+scale(VPd)+(1|fire),
                  family = binomial(link = "logit"),
                  data = FDmod)
summary(m2) ## AIC 19974.0

plot(residuals(m2),
     ylab = "residuals",
     pch = 16,
     cex = 0.75,
     col = rgb(0,0,0,0.75),
     las = 1,
     xlab = "index") ## residuals still look whack
abline(h = 0, col = "red", lwd = 2, lty = 2)
qqnorm(residuals(m2), pch = 16, cex = .75, col = rgb(0,0,0,0.75))
qqline(residuals(m2), col = "red", lwd = 2, lty = 2) ## somehow worse?

rm(m1);rm(m2);rm(pseudo.R.squared)
rm(FDmod)
## thus concludes the limited attempt to construct linear models

#### Calculating spatial autocorrelation eigen vectors for random forest ####
## models will need to be run twice to reproduce results
## once with all data, a second time with the removal of 4 largest fires (supplement)

## uncomment this code for supplemental analysis of fires minus 4 largest
#
#
# vec <- c("cameronpeak","mullen","easttroublesome","hermitspeak") ## 4 largest fires
# FireData <- FireData[!FireData$fire %in% vec,] ## subsetting out megafires
#
#
#

## spatially thinning my data to keep 1 sample per 100 x 100 m cell
dat <- FireData
r <- rast("./Predictor Rasters/Land Cover/Aspen_Binary.tif")
blank <- rast(ext(r), resolution=100, vals=NA) ## gonna expand this
dat.sp <- vect(dat, geom = c("x","y"),crs = crs(blank))
dat.cell <- extract(blank, dat.sp, cell = TRUE)
dat$cell <- dat.cell$cell
dat <- dat %>% group_by(cell) %>% sample_n(size=1) # sample one point per 100 x 100 m cell
dat$cell <- NULL

## looking at spatial autocorrelation in the data
dat$stat <- as.integer(ifelse(dat$stat == "EF", 0,1)) ## turning into integer for correlog

optimal_increment <- max(dist(dat[, c("x", "y")])) / 15 ## trying to shoot for 15 - 20 bins
## correlogram can take a few minutes to run
correlogram <- ncf::correlog(x = dat$x, y = dat$y, z = dat$stat,
                            increment = optimal_increment, ## optimal distance
                            resamp = 100, ## number of bootstrap resamples 
                            latlon = FALSE, 
                            na.rm = TRUE, 
                            quiet = FALSE)
plot(correlogram)
dist_threshold <- as.numeric(correlogram$x.intercept) 
## 33067.81 for whole dataset if the correlogram fn isn't working. This is in m
## 220238.2 for subset without large fires. Also in m

## PCNM method for spatial autocorrelation
dmat <- as.matrix(dist(cbind(dat$y, dat$x))) ## turning the coordinates of each plot into a distance matrix
pcnm_results <- pcnm(dmat, threshold = dist_threshold)
plot(pcnm_results$values)
slopes <- NA ## looking at slopes to determine when curve flattens
for(i in 1:200){
  slopes[i] <- (pcnm_results$values[i+1] - pcnm_results$values[i])/((i+1)-i) 
}
plot(slopes)
points(slopes[(1:30)], col = "red") ## 30 is sufficient

num_eigenvectors <- 30
eigen_res <- eigs_sym(as.matrix(dmat), k = num_eigenvectors)
print(eigen_res$vectors)

ncol(dat)
ncol(dat) + 30
dat[,30:59]<-eigen_res$vectors[,1:30]
colnames(dat)[30:59]

#### Random Forest ####
dat$stat <- as.factor(dat$stat)
colnames(dat)
dat_sub <- dat[c(2,9:23,28,30:59)]

level1 <- dat_sub[dat_sub$stat == levels(dat_sub$stat)[1], ]
level2 <- dat_sub[dat_sub$stat == levels(dat_sub$stat)[2], ]
set.seed(1)
EF_sample <- level1[sample(nrow(level1), 1000, replace = TRUE), ]
set.seed(1)
EH_sample <- level2[sample(nrow(level2), 1000, replace = TRUE), ]
dat_sub_training <- rbind(EF_sample, EH_sample)

train_control <- trainControl(method = "cv", number = 5, p = 0.25)
set.seed(1)
train_index <- createDataPartition(y = dat_sub_training$stat, p = 0.8, list = FALSE) ## prev. 0.75
training_set <- dat_sub_training[train_index,]
testing_set <- dat_sub_training[-train_index,]

# Tune mtry and nodesize
training_set <- as.data.frame(training_set)
tuned_model <- tune.rfsrc(stat ~ ., data = training_set,
                          nodesizeTry = c(1:9, seq(10, 100, by = 5)),
                          ntree.try = 500)
# View optimal parameters
print(tuned_model)


1-length(which(training_set$stat == 0))/length(training_set$stat) ## unbalanced withouth megafires - most lines hold

set.seed(1)
model <- train(stat~.,
               data = training_set,
               method = "rf",
               ntree = 500,
               maxnodes = 100,
               keep.inbag = TRUE,
               keep.forest = TRUE,
               maximize = TRUE,
               trControl = train_control,
               metric = "Accuracy")

model

y_hats <- predict(object = model,
                  newdata = testing_set[, -1])
## Print the accuracy
mean(as.numeric(y_hats) - as.numeric(testing_set$stat)) 
# Average difference between the two estimates
# proportional difference of < 1% -> the model did pretty good
# Sub-sampling the data to be more balanced leads to over prediction of lines failing

n <- 100

VPd.x <- matrix(data = NA, nrow = n, ncol = 51)
VPd.y <- matrix(data = NA, nrow = n, ncol = 51)

growth.x <- matrix(data = NA, nrow = n, ncol = 51)
growth.y <- matrix(data = NA, nrow = n, ncol = 51)

elev.x <- matrix(data = NA, nrow = n, ncol = 51)
elev.y <- matrix(data = NA, nrow = n, ncol = 51)

slope.x <- matrix(data = NA, nrow = n, ncol = 51)
slope.y <- matrix(data = NA, nrow = n, ncol = 51)

tpi.x <- matrix(data = NA, nrow = n, ncol = 51)
tpi.y <- matrix(data = NA, nrow = n, ncol = 51)

prop.sf.x <- matrix(data = NA, nrow = n, ncol = 51)
prop.sf.y <- matrix(data = NA, nrow = n, ncol = 51)

prop.pj.x <- matrix(data = NA, nrow = n, ncol = 51)
prop.pj.y <- matrix(data = NA, nrow = n, ncol = 51)

prop.pico.x <- matrix(data = NA, nrow = n, ncol = 51)
prop.pico.y <- matrix(data = NA, nrow = n, ncol = 51)

prop.pipo.x <- matrix(data = NA, nrow = n, ncol = 51)
prop.pipo.y <- matrix(data = NA, nrow = n, ncol = 51)

prop.psme.x <- matrix(data = NA, nrow = n, ncol = 51)
prop.psme.y <- matrix(data = NA, nrow = n, ncol = 51)

prop.oak.x <- matrix(data = NA, nrow = n, ncol = 51)
prop.oak.y <- matrix(data = NA, nrow = n, ncol = 51)

prop.shr.x <- matrix(data = NA, nrow = n, ncol = 51)
prop.shr.y <- matrix(data = NA, nrow = n, ncol = 51)

prop.oth.x <- matrix(data = NA, nrow = n, ncol = 51)
prop.oth.y <- matrix(data = NA, nrow = n, ncol = 51)

prop.asp.x <- matrix(data = NA, nrow = n, ncol = 51)
prop.asp.y <- matrix(data = NA, nrow = n, ncol = 51)

prop.grs.x <- matrix(data = NA, nrow = n, ncol = 51)
prop.grs.y <- matrix(data = NA, nrow = n, ncol = 51)

perim_m.x <- matrix(data = NA, nrow = n, ncol = 51)
perim_m.y <- matrix(data = NA, nrow = n, ncol = 51)

testing_list <- vector("list", n) 

rf.res <- matrix(data = NA, nrow = n, ncol = nrow(training_set))

y_hats <- matrix(data = NA, nrow = n, ncol = nrow(testing_set))
y_hats.diff <- NA
varImp.summary <- matrix(data = NA, nrow = length(training_set)-1, ncol = n)
varImp.names <- matrix(data = NA, nrow = length(training_set)-1, ncol = n)

balance <- NA
error <- matrix(data = NA, nrow = n, ncol = 500) ## ncol = ntree
AUC.val <- NA

## for loop for the random forest and summary data
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK")
doParallel::registerDoParallel(cl = local.cluster) ## parallel processing to run faster

for(i in 1:n){
  set.seed(i)
  EF_sample <- level1[sample(nrow(level1), 1000, replace = FALSE), ]
  set.seed(i)
  EH_sample <- level2[sample(nrow(level2), 1000, replace = FALSE), ]
  dat_sub <- rbind(EF_sample, EH_sample)
  train_index <- createDataPartition(y = dat_sub_training$stat, p = 0.8, list = FALSE) ## prev. 0.75
  training_set <- as.data.frame(dat_sub[train_index,])
  testing_set <- as.data.frame(dat_sub[-train_index,])
  balance[i] <- 1-length(which(training_set$stat == 0))/length(training_set$stat)
  
  set.seed(i)
  rf <- randomForest(stat~.,
                     data = training_set,
                     ntree = 500,
                     maxnodes = 100,
                     maximize = TRUE,
                     trControl = train_control,
                     importance = TRUE,
                     keep.forest = TRUE,
                     keep.inbag = TRUE) ## making the rf object
  y_hats[i,1:nrow(testing_set)] <- predict(object = rf, newdata = testing_set[, -1],type = "prob")[,2] ## predicted probability of holding
  testing_list[[i]] <- testing_set
  y_hats.diff[i] <- mean(round(as.numeric(y_hats[i,1:nrow(testing_set)]),0) - (as.numeric(testing_set$stat)-1))
  varImp.summary[,i] <- rf$importance[,3]
  varImp.names[,i] <- rownames(rf$importance)
  rf.res[i,] <- as.integer(rf$predicted)
  error[i,] <- rf$err.rate[,1]
  rf.roc <- suppressMessages(roc(training_set$stat, rf$votes[,2]))
  AUC.val[i] <- as.numeric(auc(rf.roc))
  
  VPd <- partialPlot(rf, training_set, x.var = VPd)
  growth <- partialPlot(rf, training_set, x.var = Growth)
  elev <- partialPlot(rf, training_set, x.var = elev)
  slope <- partialPlot(rf, training_set, x.var = slope)
  tpi <- partialPlot(rf, training_set, x.var = tpi)
  prop.sf <- partialPlot(rf, training_set, x.var = sf_prop)
  prop.pj <- partialPlot(rf, training_set, x.var = pj_prop)
  prop.pico <- partialPlot(rf, training_set, x.var = lodgepole_prop)
  prop.pipo <- partialPlot(rf, training_set, x.var = pondo_prop)
  prop.psme <- partialPlot(rf, training_set, x.var = dougfir_prop)
  prop.oak <- partialPlot(rf, training_set, x.var = gambel_prop)
  prop.shr <- partialPlot(rf, training_set, x.var = shrub_prop)
  prop.asp <- partialPlot(rf, training_set, x.var = aspen_prop)
  prop.grs <- partialPlot(rf, training_set, x.var = grass_prop)
  prop.oth <- partialPlot(rf, training_set, x.var = other_prop)
  perim_m <- partialPlot(rf, training_set, x.var = perim_m)
  
  VPd.x[i,1:length(VPd$x)] <- VPd$x
  VPd.y[i,1:length(VPd$y)] <- VPd$y
  growth.x[i,1:length(growth$x)] <- growth$x
  growth.y[i,1:length(growth$y)] <- growth$y
  elev.x[i,1:length(elev$x)] <- elev$x
  elev.y[i,1:length(elev$y)] <- elev$y
  slope.x[i,1:length(slope$x)] <- slope$x
  slope.y[i,1:length(slope$y)] <- slope$y
  tpi.x[i,1:length(tpi$x)] <- tpi$x
  tpi.y[i,1:length(tpi$y)] <- tpi$y
  prop.sf.x[i,1:length(prop.sf$x)] <- prop.sf$x
  prop.sf.y[i,1:length(prop.sf$y)] <- prop.sf$y
  prop.pj.x[i,1:length(prop.pj$x)] <- prop.pj$x
  prop.pj.y[i,1:length(prop.pj$y)] <- prop.pj$y
  prop.pipo.x[i,1:length(prop.pipo$x)] <- prop.pipo$x
  prop.pipo.y[i,1:length(prop.pipo$y)] <- prop.pipo$y
  prop.pico.x[i,1:length(prop.pico$x)] <- prop.pico$x
  prop.pico.y[i,1:length(prop.pico$y)] <- prop.pico$y
  prop.psme.x[i,1:length(prop.psme$x)] <- prop.psme$x
  prop.psme.y[i,1:length(prop.psme$y)] <- prop.psme$y
  prop.oak.x[i,1:length(prop.oak$x)] <- prop.oak$x
  prop.oak.y[i,1:length(prop.oak$y)] <- prop.oak$y
  prop.shr.x[i,1:length(prop.shr$x)] <- prop.shr$x
  prop.shr.y[i,1:length(prop.shr$y)] <- prop.shr$y
  prop.asp.x[i,1:length(prop.asp$x)] <- prop.asp$x
  prop.asp.y[i,1:length(prop.asp$y)] <- prop.asp$y
  prop.grs.x[i,1:length(prop.grs$x)] <- prop.grs$x
  prop.grs.y[i,1:length(prop.grs$y)] <- prop.grs$y
  prop.oth.x[i,1:length(prop.oth$x)] <- prop.oth$x
  prop.oth.y[i,1:length(prop.oth$y)] <- prop.oth$y
  perim_m.x[i,1:length(perim_m$x)] <- perim_m$x
  perim_m.y[i,1:length(perim_m$y)] <- perim_m$y
  gc()
  
  progress <- i/n*100
  if (progress %% 5 == 0) {
    print(paste(progress, "% done", sep = ""))
  }
}
parallel::stopCluster(cl = local.cluster)

## pred vs obs plots
y_hats.diff <- y_hats.diff*100 ## converting to %
plot(x = 1:length(y_hats.diff), y = y_hats.diff,
     pch = 16,
     xlab = "model run",
     ylim = c(min(y_hats.diff)-10,max(y_hats.diff)+10),
     las = 1,
     ylab = "% Difference in Predicted vs.Observed",
     cex = 1)
abline(h = mean(y_hats.diff), col="firebrick4", lty = 2)
text(x = 30, y = 10, paste("Average difference = ", round(mean(y_hats.diff), digits = 1),"%", sep = "")) 

mean(balance);min(balance);max(balance)
# 0.5 - 0.9232063
# 0.5 - 0.8952381
# 0.5 - 0.9507937 far less balance on the no mega fire subset

error.mean <- apply(error,2,mean, na.rm = TRUE)
min(error, na.rm = TRUE);max(error, na.rm = TRUE)
plot(error.mean, type = "n",
     ylim = c(0,max(error, na.rm = TRUE)+0.1),
     xlab = "Tree",
     ylab = "Error")
for(i in 1:100){
  lines(error[i,], col = rgb(0,0,0,alpha = 0.25))
}
lines(error.mean, type = "l", col = "firebrick", lty = 2, lwd= 2)
error.mean[500]*100
text(x = 300, y = 0.3, paste("Average error = ", round(error.mean[500]*100, digits = 1),"%", sep = "")) 

mean(AUC.val);min(AUC.val);max(AUC.val)
# 0.866342 0.8435496
# 0.8501945 0.7550632
# 0.8845094 0.9114035

## VarImp Plot
varImp.names[c(17:46),1]
varImp.plotting <- data.frame(name = c(varImp.names[c(1:16),1],"spatial"),
                              mean = c(apply(varImp.summary[c(1:16),],1,mean,na.rm = TRUE),mean(varImp.summary[c(17:46),],na.rm = TRUE)),
                              min = c(apply(varImp.summary[c(1:16),],1,min,na.rm = TRUE),min(varImp.summary[c(17:46),],na.rm = TRUE)),
                              max = c(apply(varImp.summary[c(1:16),],1,max,na.rm = TRUE),max( varImp.summary[c(17:46),],na.rm = TRUE)))
# varImp.plotting$name <- c("VPd","Elevation","Slope", "TPI","Aspen","Spruce/Fir",
#                           "Pinyon/Juniper","Lodgepole","Ponderosa","Douglas/Fir",
#                           "Oak","Shrubland","Grassland","Other","Growth", "Spatial")

varImp.plotting <- varImp.plotting[order(varImp.plotting$mean, decreasing = FALSE),]

len <- as.integer(nrow(varImp.plotting))
par(mfrow = c(1,1), oma = c(0,3,0,0))
plot(varImp.plotting$mean,
     ylim = c(0,len),
     xlim = c(0,max(varImp.plotting$max,na.rm = TRUE)),
     las = 1,
     type = "n",
     ylab = "",
     yaxt = "n",
     xlab = "Variable Importance")
axis(2, at = c(1:len), labels = varImp.plotting$name, cex.axis = 1, las = 2)
points(x = varImp.plotting$mean,y = 1:len, col = "black", cex = 1, pch = 16)
segments(x0 = varImp.plotting$min, y0 = 1:len, x1 = varImp.plotting$max, y1 = 1:len, col = "black", lwd = 1.5)

## partial dependence plots
FD <- dat_sub
FD$LineInt <- as.integer(FD$stat)-1
FD$LineInt[FD$LineInt == 0] <- -0.25
FD$LineInt[FD$LineInt == 1] <- 1.25

# VPd.y <- 1-(1/(1+exp(-VPd.y)))*2
# growth.y <- 1-(1/(1+exp(-growth.y)))*2
# elev.y <- 1-(1/(1+exp(-elev.y)))*2
# slope.y <- 1-(1/(1+exp(-slope.y)))*2
# tpi.y <- 1-(1/(1+exp(-tpi.y)))*2
# perim_m.y <- 1-(1/(1+exp(-perim_m.y)))*2
# prop.sf.y <- 1-(1/(1+exp(-prop.sf.y)))*2
# prop.pj.y <- 1-(1/(1+exp(-prop.pj.y)))*2
# prop.pico.y <- 1-(1/(1+exp(-prop.pico.y)))*2
# prop.pipo.y <- 1-(1/(1+exp(-prop.pipo.y)))*2
# prop.psme.y <- 1-(1/(1+exp(-prop.psme.y)))*2
# prop.oak.y <- 1-(1/(1+exp(-prop.oak.y)))*2
# prop.shr.y <- 1-(1/(1+exp(-prop.shr.y)))*2
# prop.oth.y <- 1-(1/(1+exp(-prop.oth.y)))*2
# prop.asp.y <- 1-(1/(1+exp(-prop.asp.y)))*2
# prop.grs.y <- 1-(1/(1+exp(-prop.grs.y)))*2
gc()

par(mfrow = c(4,4))

plot(growth.x[1,], growth.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Fire Growth (ha)",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(growth.x[i,], growth.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$Growth, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
growth.x.mean <- apply(growth.x,2,mean, na.rm = T)
growth.y.mean <- apply(growth.y,2,mean,, na.rm = T)
lo <- loess(growth.y.mean~growth.x.mean, na.rm = TRUE)
lines(y = predict(lo), x = growth.x.mean[1:length(predict(lo))], col = "red", lwd = 2)
# fire growth over 5947.384 associated with failure

plot(VPd.x[1,], VPd.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "VPd Difference (kPa)",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(VPd.x[i,], VPd.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$VPd, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
VPd.x.mean <- apply(VPd.x,2,mean)
VPd.y.mean <- apply(VPd.y,2,mean)
lo <- loess(VPd.y.mean~VPd.x.mean)
lines(y = predict(lo), x = VPd.x.mean, col = "red", lwd = 2)

plot(elev.x[1,], elev.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Elevation (m)",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(elev.x[i,], elev.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$elev, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
elev.x.mean <- apply(elev.x,2,mean)
elev.y.mean <- apply(elev.y,2,mean)
lo <- loess(elev.y.mean~elev.x.mean)
lines(y = predict(lo), x = elev.x.mean, col = "red", lwd = 2)

plot(slope.x[1,], slope.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Slope (degrees)",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(slope.x[i,], slope.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$slope, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
slope.x.mean <- apply(slope.x,2,mean)
slope.y.mean <- apply(slope.y,2,mean)
lo <- loess(slope.y.mean~slope.x.mean)
lines(y = predict(lo), x = slope.x.mean, col = "red", lwd = 2)

plot(tpi.x[1,], tpi.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Topographic Position Index",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(tpi.x[i,], tpi.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$tpi, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
tpi.x.mean <- apply(tpi.x,2,mean)
tpi.y.mean <- apply(tpi.y,2,mean)
lo <- loess(tpi.y.mean~tpi.x.mean)
lines(y = predict(lo), x = tpi.x.mean, col = "red", lwd = 2)

plot(perim_m.x[1,], perim_m.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Length of Fire Line (m)",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(perim_m.x[i,], perim_m.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$perim_m, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
perim_m.x.mean <- apply(perim_m.x,2,mean)
perim_m.y.mean <- apply(perim_m.y,2,mean)
lo <- loess(perim_m.y.mean~perim_m.x.mean)
lines(y = predict(lo), x = perim_m.x.mean, col = "red", lwd = 2)

## community types
plot(prop.sf.x[1,], prop.sf.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Proportion Spruce/Fir",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(prop.sf.x[i,], prop.sf.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$sf_prop, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
prop.sf.x.mean <- apply(prop.sf.x,2,mean)
prop.sf.y.mean <- apply(prop.sf.y,2,mean)
lo <- loess(prop.sf.y.mean~prop.sf.x.mean)
lines(y = predict(lo), x = prop.sf.x.mean, col = "red", lwd = 2)

plot(prop.pico.x[1,], prop.pico.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Proportion Lodgepole",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(prop.pico.x[i,], prop.pico.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$lodgepole_prop, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
prop.pico.x.mean <- apply(prop.pico.x,2,mean)
prop.pico.y.mean <- apply(prop.pico.y,2,mean)
lo <- loess(prop.pico.y.mean~prop.pico.x.mean)
lines(y = predict(lo), x = prop.pico.x.mean, col = "red", lwd = 2)

plot(prop.asp.x[1,], prop.asp.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Proportion Aspen",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(prop.asp.x[i,], prop.asp.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$aspen_prop, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
prop.asp.x.mean <- apply(prop.asp.x,2,mean)
prop.asp.y.mean <- apply(prop.asp.y,2,mean)
lo <- loess(prop.asp.y.mean~prop.asp.x.mean)
lines(y = predict(lo), x = prop.asp.x.mean, col = "red", lwd = 2)

plot(prop.psme.x[1,], prop.psme.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Proportion Douglas Fir",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(prop.psme.x[i,], prop.psme.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$dougfir_prop, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
prop.psme.x.mean <- apply(prop.psme.x,2,mean)
prop.psme.y.mean <- apply(prop.psme.y,2,mean)
lo <- loess(prop.psme.y.mean~prop.psme.x.mean)
lines(y = predict(lo), x = prop.psme.x.mean, col = "red", lwd = 2)

plot(prop.pipo.x[1,], prop.pipo.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Proportion Ponderosa",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(prop.pipo.x[i,], prop.pipo.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$pondo_prop, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
prop.pipo.x.mean <- apply(prop.pipo.x,2,mean)
prop.pipo.y.mean <- apply(prop.pipo.y,2,mean)
lo <- loess(prop.pipo.y.mean~prop.pipo.x.mean)
lines(y = predict(lo), x = prop.pipo.x.mean, col = "red", lwd = 2)

plot(prop.pj.x[1,], prop.pj.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Proportion Pinyon/Juniper",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(prop.pj.x[i,], prop.pj.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$pj_prop, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
prop.pj.x.mean <- apply(prop.pj.x,2,mean)
prop.pj.y.mean <- apply(prop.pj.y,2,mean)
lo <- loess(prop.pj.y.mean~prop.pj.x.mean)
lines(y = predict(lo), x = prop.pj.x.mean, col = "red", lwd = 2)

plot(prop.oak.x[1,], prop.oak.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Proportion Oak",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(prop.oak.x[i,], prop.oak.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$gambel_prop, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
prop.oak.x.mean <- apply(prop.oak.x,2,mean)
prop.oak.y.mean <- apply(prop.oak.y,2,mean)
lo <- loess(prop.oak.y.mean~prop.oak.x.mean)
lines(y = predict(lo), x = prop.oak.x.mean, col = "red", lwd = 2)

plot(prop.shr.x[1,], prop.shr.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Proportion Shrubland",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(prop.shr.x[i,], prop.shr.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$shrub_prop, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
prop.shr.x.mean <- apply(prop.shr.x,2,mean)
prop.shr.y.mean <- apply(prop.shr.y,2,mean)
lo <- loess(prop.shr.y.mean~prop.shr.x.mean)
lines(y = predict(lo), x = prop.shr.x.mean, col = "red", lwd = 2)

plot(prop.grs.x[1,], prop.grs.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Proportion Grasslands",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(prop.grs.x[i,], prop.grs.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$grass_prop, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
prop.grs.x.mean <- apply(prop.grs.x,2,mean)
prop.grs.y.mean <- apply(prop.grs.y,2,mean)
lo <- loess(prop.grs.y.mean~prop.grs.x.mean, na.rm = T)
lines(y = predict(lo), x = prop.grs.x.mean[1:length(predict(lo))], col = "red", lwd = 2)

plot(prop.oth.x[1,], prop.oth.y[1,],
     type = "l",
     ylim = c(-0.25,1.25),
     col = rgb(0,0,0,0.25),
     main = "",
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 1,
     xlab = "Proportion Other",
     # ylab = "(Logit of Probability of Line Holding)/2")
     ylab = "Line Status")
for(i in 2:50)(
  lines(prop.oth.x[i,], prop.oth.y[i,], col = rgb(0,0,0,0.25))
)
abline(h = 0.5, lty = 2, lwd = 2)
points(x = FD$other_prop, y = jitter(FD$LineInt, factor = .25),
       cex = 0.5,
       col = rgb(0,0,0,0.25),
       pch = 16)
prop.oth.x.mean <- apply(prop.oth.x,2,mean)
prop.oth.y.mean <- apply(prop.oth.y,2,mean)
lo <- loess(prop.oth.y.mean~prop.oth.x.mean)
lines(y = predict(lo), x = prop.oth.x.mean, col = "red", lwd = 2)

## linear model of predicted probability of line holding
par(mfrow = c(1,1))
## creating empty dataframes/vectors
p.val <- data.frame(matrix(ncol = 17, nrow = 100))
colnames(p.val) <- c("intercept", colnames(testing_list[[1]])[2:17])
est <- data.frame(matrix(ncol = 17, nrow = 100))
colnames(est) <- c("intercept", colnames(testing_list[[1]])[2:17])
r.sq <- NA

newdata_list <- vector("list", n) ## creating list to store 'newdata' holds all other predictors at average while varying partial effects of predictor of interest
preds_list <- vector("list", n) ## output list for predictors
vec <- seq(1:16) ## number of predictors
## this for loop goes through each of i:100 training_set dfs to make linear models
## after models are constructed the partial effects of j:16 predictors are pulled out and stored

for(i in 1:n){
  mod <- lm(y_hats[i,] ~ ., data = testing_list[[i]][2:17])
  mod.sum <- summary(mod)
  p.val[i,] <- mod.sum$coefficients[,4]
  est[i,] <- mod.sum$coefficients[,1]
  r.sq[i] <- mod.sum$r.squared
  
  for(j in 1:length(vec)){
    newdata <- data.frame(matrix(ncol = 16, nrow = 30))
    colnames(newdata) <- c(colnames(testing_list[[i]])[2:17])
    newdata[,j] <- seq(min(testing_list[[i]][j+1]),max(testing_list[[i]][j+1]), length.out = 30)
    columnMeans <- colMeans(testing_list[[i]][2:17])
    newdata[1:30,vec[-j]] <- rep(columnMeans[vec[-j]],each = 30)
    newdata_list[[i]][[j]] <- newdata
    preds <- predict(mod, newdata, type="response")
    preds_list[[i]][[j]] <- preds
  }
}

hist(r.sq) ## looking at the distribution of r.sq
for(i in 1:ncol(p.val)){
  ifelse(length(which(p.val[,i] < 0.05))>5, print(colnames(p.val)[i]), NA)
} ## if p values < 0.05 for more than 5 of 100 model runs, there is likely a difference in the distributions 

newdata_avg <- data.frame(matrix(ncol = 30, nrow = 100)) ## storing newdata for each of j:16 predictors
preds_avg <- data.frame(matrix(ncol = 30, nrow = 100)) ## storing predictions as well

par(mfrow = c(4,4)) ## reset the plotting to 4 x 4
## for loop to plot each of the data and lines for the partial effect of each predictor variable (j)
for(j in 1:length(vec)){
  plot(y_hats[1,] ~ as.numeric(testing_list[[1]][,j+1]),
       ylab = "predicted line status",
       las = 1,
       xlab = colnames(testing_list[[1]])[j+1],
       pch = 16,
       col = adjustcolor("black",alpha.f = 0.01))
  for(i in 1:n){
    points(y_hats[i,] ~ as.numeric(testing_list[[i]][,j+1]),
           pch = 16,
           col = adjustcolor("black",alpha.f = 0.01))
    lines(newdata_list[[i]][[j]][,j], preds_list[[i]][[j]], lty=1, col = adjustcolor("black",alpha.f = 0.1))
    newdata_avg[i,] <- newdata_list[[i]][[j]][,j]
    preds_avg[i,] <- preds_list[[i]][[j]]
  }
  abline(h = 0.5, lty = 2)
  lines(apply(newdata_avg,2,mean), apply(preds_avg,2,mean), col = "red")
}

rm(list = ls()) ## cleaning global env
gc()

#### Nonparametric Models ####
SR <- vect("./Boundary/SouthernRockyBoundary_10kmBuff.shp") ## loading in the study area boundary, may or may not need to re:load based on where in script you are starting
SR_Fires <- vect("./NIFC Polygons/SR_Fires.shp")
SR_FLs <- vect("./NIFC Lines/SR_FLs.shp")
Fires_add60 <- buffer(SR_Fires, 60)
FLs_Engaged <- terra::intersect(SR_FLs, Fires_add60)

## pulling in the landcover rasters again to count the number of pixels
temp <- list.files(path = "./Predictor Rasters/Land Cover/",pattern="*.tif") ## creating a vector that has all the files in the working directory with .tif extensions
for(i in 1:length(temp)) {
  path <- paste("./Predictor Rasters/Land Cover/", temp[i], sep = "") ## specifying the relative pathway for assign
  assign(temp[i], terra::rast(path)) ## assigning the rasters
} ## loading in the raters I want

sr_masked <- mask(Aspen_Binary.tif, SR) ## masking the layer to the SR
tots <- as.integer(freq(sr_masked, value = 1)[3]) ## recording the number of cells == 1
fires <- mask(sr_masked, SR_Fires) ## masking furhter to the burned areas
burns <- as.integer(freq(fires, value = 1)[3]) ## recording the number of cells == 1
FLs <- mask(sr_masked, FLs_Engaged, touches = TRUE) ## final mask to the fire lines
Fls <- as.integer(freq(FLs, value = 1)[3]) ## recording these cells
gc()

## repeating, but making sure that further additions are binded to original outputs
## dougfir
sr_masked <- mask(DougFir_Binary.tif, SR)
tots.added <- as.integer(freq(sr_masked, value = 1)[3])
tots <- c(tots, tots.added)
fires <- mask(sr_masked, SR_Fires)
burns.added <- as.integer(freq(fires, value = 1)[3])
burns <- c(burns, burns.added)
FLs <- mask(sr_masked, FLs_Engaged, touches = TRUE)
Fls.added <- as.integer(freq(FLs, value = 1)[3]) 
Fls <- c(Fls, Fls.added)
gc()

## gambel
sr_masked <- mask(Gambel_Binary.tif, SR)
tots.added <- as.integer(freq(sr_masked, value = 1)[3])
tots <- c(tots, tots.added)
fires <- mask(sr_masked, SR_Fires)
burns.added <- as.integer(freq(fires, value = 1)[3])
burns <- c(burns, burns.added)
FLs <- mask(sr_masked, FLs_Engaged, touches = TRUE)
Fls.added <- as.integer(freq(FLs, value = 1)[3])
Fls <- c(Fls, Fls.added)
gc()

## grass
sr_masked <- mask(Grass_Binary.tif, SR)
tots.added <- as.integer(freq(sr_masked, value = 1)[3])
tots <- c(tots, tots.added)
fires <- mask(sr_masked, SR_Fires)
burns.added <- as.integer(freq(fires, value = 1)[3])
burns <- c(burns, burns.added)
FLs <- mask(sr_masked, FLs_Engaged, touches = TRUE)
Fls.added <- as.integer(freq(FLs, value = 1)[3])
Fls <- c(Fls, Fls.added)
gc()

## lodgepole
sr_masked <- mask(Lodgepole_Binary.tif, SR)
tots.added <- as.integer(freq(sr_masked, value = 1)[3])
tots <- c(tots, tots.added)
fires <- mask(sr_masked, SR_Fires)
burns.added <- as.integer(freq(fires, value = 1)[3])
burns <- c(burns, burns.added)
FLs <- mask(sr_masked, FLs_Engaged, touches = TRUE)
Fls.added <- as.integer(freq(FLs, value = 1)[3])
Fls <- c(Fls, Fls.added)
gc()

## other
sr_masked <- mask(Other_Binary.tif, SR)
tots.added <- as.integer(freq(sr_masked, value = 1)[3])
tots <- c(tots, tots.added)
fires <- mask(sr_masked, SR_Fires)
burns.added <- as.integer(freq(fires, value = 1)[3])
burns <- c(burns, burns.added)
FLs <- mask(sr_masked, FLs_Engaged, touches = TRUE)
Fls.added <- as.integer(freq(FLs, value = 1)[3])
Fls <- c(Fls, Fls.added)
gc()

## PJ
sr_masked <- mask(PJ_Binary.tif, SR)
tots.added <- as.integer(freq(sr_masked, value = 1)[3])
tots <- c(tots, tots.added)
fires <- mask(sr_masked, SR_Fires)
burns.added <- as.integer(freq(fires, value = 1)[3])
burns <- c(burns, burns.added)
FLs <- mask(sr_masked, FLs_Engaged, touches = TRUE)
Fls.added <- as.integer(freq(FLs, value = 1)[3])
Fls <- c(Fls, Fls.added)
gc()

## ponderosa
sr_masked <- mask(Ponderosa_Binary.tif, SR)
tots.added <- as.integer(freq(sr_masked, value = 1)[3])
tots <- c(tots, tots.added)
fires <- mask(sr_masked, SR_Fires)
burns.added <- as.integer(freq(fires, value = 1)[3])
burns <- c(burns, burns.added)
FLs <- mask(sr_masked, FLs_Engaged, touches = TRUE)
Fls.added <- as.integer(freq(FLs, value = 1)[3]) 
Fls <- c(Fls, Fls.added)
gc()

## SF
sr_masked <- mask(SF_Binary.tif, SR)
tots.added <- as.integer(freq(sr_masked, value = 1)[3])
tots <- c(tots, tots.added)
fires <- mask(sr_masked, SR_Fires)
burns.added <- as.integer(freq(fires, value = 1)[3])
burns <- c(burns, burns.added)
FLs <- mask(sr_masked, FLs_Engaged, touches = TRUE)
Fls.added <- as.integer(freq(FLs, value = 1)[3])
Fls <- c(Fls, Fls.added)
gc()

## shrub
sr_masked <- mask(Shrub_Binary.tif, SR)
tots.added <- as.integer(freq(sr_masked, value = 1)[3])
tots <- c(tots, tots.added)
fires <- mask(sr_masked, SR_Fires)
burns.added <- as.integer(freq(fires, value = 1)[3])
burns <- c(burns, burns.added)
FLs <- mask(sr_masked, FLs_Engaged, touches = TRUE)
Fls.added <- as.integer(freq(FLs, value = 1)[3]) 
Fls <- c(Fls, Fls.added)
gc()


ID <- tolower(gsub("_Binary.tif", x = temp, ""))
cellcounts <- data.frame(name = ID,
                         tot = tots,
                         burn = burns,
                         Fls = Fls)
write.csv(cellcounts, "./Results/cellcounts.csv")
rm(list = temp)

## statistical analysis
cellcounts <- read.csv("./Results/cellcounts.csv")

## burned vs landscape
burned <- cellcounts$burn
unburned <- cellcounts$tot - cellcounts$burn
BurnedMat <- matrix(c(burned, unburned), nrow = 10)
colnames(BurnedMat) <- c("burned", "unburned")
rownames(BurnedMat) <- cellcounts$name
chisq.test(BurnedMat) ## difference in representation in burns

Burned.Res <- data.frame(name = cellcounts$name[-1],
                         est = NA,
                         lwr = NA,
                         upr = NA)
for(i in 2:nrow(BurnedMat)){
  mod <- fisher.test(BurnedMat[c(i,1),])
  Burned.Res$est[i-1] <- mod$estimate
  Burned.Res$lwr[i-1] <- mod$conf.int[1]
  Burned.Res$upr[i-1] <- mod$conf.int[2]
}

## fire lines vs burned
firelines <- cellcounts$Fls
nolines <- cellcounts$burn - cellcounts$Fls
LinesMat <- matrix(c(firelines, nolines), nrow = 10)
colnames(LinesMat) <- c("lines", "nolines")
rownames(LinesMat) <- cellcounts$name
chisq.test(LinesMat) ## difference in representation in burns

Lines.Res <- data.frame(name = cellcounts$name[-1],
                        est = NA,
                        lwr = NA,
                        upr = NA)
for(i in 2:nrow(LinesMat)){
  mod <- fisher.test(LinesMat[c(i,1),])
  Lines.Res$est[i-1] <- mod$estimate
  Lines.Res$lwr[i-1] <- mod$conf.int[1]
  Lines.Res$upr[i-1] <- mod$conf.int[2]
}

NonParametricResults <- cbind(Burned.Res,Lines.Res)
NonParametricResults <- NonParametricResults[,-5]
colnames(NonParametricResults) <- c("name", "burn.est", "burn.lwr", "burn.upr","fl.est","fl.lwr","fl.upr")
write.csv(NonParametricResults,"./Results/NonParametricResults.csv")

## Plotting Non Parametric Models
par(mfrow = c(1,1),oma = c(0,2.5,0,0))

max(Burned.Res$upr);max(Lines.Res$upr)

plot(x = Burned.Res$est,
     y = c(1:9),
     xlim = c(0,6),
     las = 1,
     ylab = "",
     type = "n",
     yaxt = "n",
     xlab = "odds")
axis(2, at = c(1:9), labels = rev(Burned.Res$name), cex.axis = 1, las = 2)

points(x = rev(Lines.Res$est),y = c(1:9), cex = 1, col = "red", pch = 16)
segments(x0 = rev(Lines.Res$lwr), y0 = c(1:9), x1 = rev(Lines.Res$upr), y1 = c(1:9), col = "red",lwd = 1.5)

points(x = rev(Burned.Res$est),y = c(1:9), cex = 1, pch = 16)
segments(x0 = rev(Burned.Res$lwr), y0 = c(1:9), x1 = rev(Burned.Res$upr), y1 = c(1:9),lwd = 1.5)
abline(v = 1, lty = 2)

legend("bottomright",
       pch = 15,
       col = c("black", "red"),
       bty = "n",
       legend = c("Burn","Fire Line"))
par(oma = c(0,0,0,0))

rm(list = ls())
gc()

#### Fire Line Status by Land Cover as a Function of Growth ####
FireData <- read.csv("./Results/CleanedFireLines.csv")

## uncomment this code for supplemental analysis of 
# vec <- c("cameronpeak","mullen","easttroublesome","hermitspeak") ## 4 largest fires
# FireData <- FireData[!FireData$fire %in% vec,] ## subsetting out megafires

par(mfrow = c(1,1))

veg_types <- gsub("_prop", x = (colnames(FireData)[13:22]),replace = "")

# Function to calculate stats for a given LineStat value
calc_stats <- function(data, stat, veg_types) {
  sapply(veg_types, function(v) {
    prop_col <- paste0(v,"_prop")
    subset <- data$Growth[data$stat == stat & data[[prop_col]] >= 0.5]
    c(mean = mean(subset), 
      se = sd(subset) / sqrt(length(subset)))
  })
}

# Calculate for failed (LineStat == 0) and held (LineStat == 1)
failed_stats <- calc_stats(FireData, "EF", veg_types)
held_stats <- calc_stats(FireData, "EH", veg_types)

# Create the data frame
gr.stat <- data.frame(
  failed = failed_stats["mean", ],
  held = held_stats["mean", ],
  SE.failed = failed_stats["se", ],
  SE.held = held_stats["se", ]
)

gr.stat$lwr.fail <- gr.stat$failed - (1.96*gr.stat$SE.failed)
gr.stat$lwr.held <- gr.stat$held - (1.96*gr.stat$SE.held)

gr.stat$upr.fail <- gr.stat$failed + (1.96*gr.stat$SE.failed)
gr.stat$upr.held <- gr.stat$held + (1.96*gr.stat$SE.held)

par(oma = c(1.5,1,0,0))
max(gr.stat)
plot(gr.stat$failed,
     ylim = c(0,max(gr.stat)+500), 
     las = 1,
     ylab = "",
     type = "n",
     xaxt = "n",
     xlab = "")
mtext("Fire Growth (ha)", side = 2, line = 3.5)
axis(1, at = c(1:10), labels = rownames(gr.stat), cex.axis = 1, las = 2)

points(x = c(0.9,1.9,2.9,3.9,4.9,5.9,6.9,7.9,8.9,9.9),y = gr.stat$failed, col = "darkgoldenrod", cex = 1, pch = 16)
segments(x0 = c(0.9,1.9,2.9,3.9,4.9,5.9,6.9,7.9,8.9,9.9), y0 = gr.stat$lwr.fail, x1 = c(0.9,1.9,2.9,3.9,4.9,5.9,6.9,7.9,8.9,9.9), y1 = gr.stat$upr.fail, col = "darkgoldenrod",lwd = 1.5)
points(x = c(1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1,9.1,10.1),y = gr.stat$held, col = "navy", cex = 1, pch = 16)
segments(x0 = c(1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1,9.1,10.1), y0 = gr.stat$lwr.held, x1 = c(1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1,9.1,10.1), y1 = gr.stat$upr.held, col = "navy",lwd = 1.5)

legend("topright",
       pch = 16,
       col = c("darkgoldenrod", "navy"),
       bty = "n",
       legend = c("Failed Lines","Held Lines"))

## t.tests
t.test(seq(1:10),seq(2:20))$p.value

calc_stats <- function(data, stat1, stat2, veg_types) {
  sapply(veg_types, function(v) {
    prop_col <- paste0(v,"_prop")
    subset1 <- data$Growth[data$stat == stat1 & data[[prop_col]] >= 0.5]
    subset2 <- data$Growth[data$stat == stat2 & data[[prop_col]] >= 0.5]
    c(estimate.x = t.test(subset1, subset2)$estimate[1],
      estimate.y = t.test(subset1, subset2)$estimate[2],
      conf.int.lwr = t.test(subset1, subset2)$conf.int[1],
      conf.int.upr = t.test(subset1, subset2)$conf.int[2],
      p.val = t.test(subset1, subset2)$p.value)
  })
}

t.test_results <- as.data.frame(calc_stats(FireData, "EF","EH", veg_types))
colnames(t.test_results)[t.test_results[5,] > 0.05]
# "aspen" "grass" "other" "shrub" 
## no significant differences in growth between failed and held lines

## no large fires subset results
# "aspen"     "dougfir"   "gambel"    "lodgepole" "pj"        "pondo"     "sf"   
## no significant differences in growth between failed and held lines
## only significant differences, are grass, other, shrub


rm(gr.stat)
