#### Begin - Loading Dependencies, Set Working Directory ####
library(terra) ## needed for all spatial analysis
library(exactextractr) ## faster extract function than terra::extract
library(dplyr) ## less error associated with rounding compared to unique()
setwd("D:/Aspen Firelines/data") ## set working directory to the data folder


#### Loading in the Fire Lines Data ####
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

#### Loading Fire Polygon Data from NIFC (2019-2023) ####
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

#### Loading Processed Spatial Data ####
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

Engaged_Lines <- do.call(rbind, Engaged_Lines_list);gc()
table(Engaged_Lines$year) ## woah, way more data here
table(Engaged_Lines$Stat) ## also here
length(unique(Engaged_Lines$IncidentNa.1)) ## more fires too
table(Engaged_Lines$IncidentNa.1) ## a few errors in the fire names

## correcting a few names
Engaged_Lines$IncidentNa.1[Engaged_Lines$IncidentNa.1 == "calfcanyon"] <- "hermitspeak"
Engaged_Lines$IncidentNa.1[Engaged_Lines$IncidentNa.1 == "middleforkfire"] <- "middlefork"
Engaged_Lines$IncidentNa.1[Engaged_Lines$IncidentNa.1 == "nidnight"] <- "midnight"
Engaged_Lines$IncidentNa.1[Engaged_Lines$IncidentNa.1 == "mullenfire"] <- "mullen"
length(unique(Engaged_Lines$IncidentNa.1)) ## 52 total fires

perim_check <- data.frame(name = SR_Fires$IncidentNa,
                          perimeter = perim(SR_Fires))
perim_check$name[perim_check$name == "calfcanyon"] <- "hermitspeak"
perim_check$name[perim_check$name == "middleforkfire"] <- "middlefork"
perim_check$name[perim_check$name == "nidnight"] <- "midnight"
perim_check$name[perim_check$name == "mullenfire"] <- "mullen"

perim_check <- perim_check[order(perim_check$name),]
perim_check <- perim_check[!duplicated(perim_check$name),] ## removing duplicated (lower area)

table(Engaged_Lines$perim_EH >= perim_check$perimeter[match(Engaged_Lines$IncidentNa.1, perim_check$name)])
## this table is telling me there are 4 observations where the line perimeter >= fire perimeter

Engaged_Lines$DeleteThis <- Engaged_Lines$perim_EH >= perim_check$perimeter[match(Engaged_Lines$IncidentNa.1, perim_check$name)]
Engaged_Lines <- Engaged_Lines[Engaged_Lines$DeleteThis == FALSE,] ## only keeping FALSE values

Engaged_Lines$DeleteThis <- Engaged_Lines$perim_EH >= (perim_check$perimeter[match(Engaged_Lines$IncidentNa.1, perim_check$name)]*.8)
table(Engaged_Lines$DeleteThis) ## removing 4 more values that the line is 0.80 the sie of the fire perimeter
Engaged_Lines <- Engaged_Lines[Engaged_Lines$DeleteThis == FALSE,] ## only keeping FALSE values

str(Engaged_Lines$CreateDate)
Engaged_Lines$CreateDate <- as.Date(Engaged_Lines$CreateDate)
plot(Engaged_Lines$CreateDate)
max(Engaged_Lines$CreateDate, na.rm = TRUE)
table(is.na(Engaged_Lines$CreateDate))
Engaged_Lines <- Engaged_Lines[complete.cases(Engaged_Lines$CreateDate),] ## removing the data we have no date for

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
                             "201911","202011","202111","202211","202311")
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
}

hist(Engaged_Lines$VPd) ## well that looks about right, soooo
## next is addding wind, even though I don't use it
## then fire growth
## then on to the actual models I guess

