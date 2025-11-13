#### Begin - Loading Dependencies, Set WorkingDirector ####
library(terra) ## needed for all spatial analysis
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

#### Loading Environmental Data ####
SR <- vect("./Boundary/SouthernRockyBoundary_10kmBuff.shp") ## loading in the study area boundary, may or may not need to re:load based on where in script you are starting

temp <- list.files(path = "./Predictor Rasters/Land Cover/",pattern="*.tif") ## creating a vector that has all the files in the working directory with .tif extensions
for(i in 1:length(temp)) {
  path <- paste("./Predictor Rasters/Land Cover/", temp[i], sep = "") ## specifying the relative pathway for assign
  assign(temp[i], terra::rast(path)) ## assigning the rasters
} ## loading in the rsaters I want

crs(DouglasFir_Binary.tif)


#### Extracting Environmental Variables for Fire Lines ####


