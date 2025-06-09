# ---
# title: CONI roost habitat selection - extracting landscape scale covariate with GEE
# author: Elly Knight
# created: November 9, 2021
# updated: June 9, 2025
# ---

# Preamble - load packages
library(tidyverse) # data wrangling
library(sf) # working with shps
library(data.table) # list handling
library(lubridate) # date manipulation
library(rgee) # connect to GEE

#1. Initialize rgee----
#ee_install()
ee_Initialize()
ee_check()

#2. Settings---
#Buffer radius
rad <- 300

#3. Write functions----

#3a. Function to add property with time in milliseconds
add_date<-function(feature) {
  date <- ee$Date(ee$String(feature$get("Date")))$millis()
  feature$set(list(date_millis=date))
}

#3b. Function to buffer points
buffer_points <- function(feature){
  properties <- c("Date", "ID", "Radius", "Type", "date_millis", "ptID", "ptIDn", "X", "Y", "timestamp", "uniq")
  new_geometry <- feature$geometry()$buffer(rad)
  ee$Feature(new_geometry)$copyProperties(feature, properties)
}

#4. Import data----
raw <- read.csv("Interim/CONIMCP_Landscape.csv", header = T)
trackingdata <- raw |> 
  st_as_sf(coords = c('X','Y'), crs = 3857) |> 
  st_transform(crs=4326) |> 
  st_coordinates() |> 
  cbind(raw |> 
          dplyr::select(-X, -Y)) |> 
  rename(timestamp=DateTime) |> 
  arrange(ptID) |> 
  fill(timestamp, .direction="down") |> 
  group_by(ptID) |> 
  mutate(ptIDn = paste0(ptID, "-", row_number())) |> 
  ungroup()
  
#5. Wrangle data----
trackingdata$Date <- as.POSIXct(trackingdata$timestamp, format = "%Y-%m-%d %H:%M:%S", tz="UTC") 
trackingdata$Date <- as.factor(trackingdata$Date)
trackingdata$year <- year(ymd_hms(trackingdata$timestamp))
trackingdata$Date <- sub(" ", "T", trackingdata$Date) #Put in a format that can be read by javascript

#Sample a few points for testing----
#trackingdata <- trackingdata |> sample_n(1000)

#6. Set up loop to go through each year----
years <- unique(trackingdata$year)

data.year <- list()
data.out <- list()
for(i in 1:length(years)){
  
  year.i <- years[i]
  
  data.i <- trackingdata |> 
    dplyr::filter(year==year.i) |> 
    mutate(row = row_number(),
           n = ceiling(row/1000))
  
  loops <- max(data.i$n)
  
  for(j in 1:loops){
    
    start_time <- Sys.time()
    
    data.j <- data.i |> 
      dplyr::filter(n==j)
    
    #7. Create sf object----
    datasf <- st_as_sf(data.j, coords = c('X','Y'), crs = 4326)

    #8. Send data to GEE----
    data <- sf_as_ee(datasf)
    
    #9. Transform day into milliseconds----
    data <-data$map(add_date)
    
    #10. Buffer points----
    data.buff <- data$map(buffer_points)

    #11. Load EVI image collection----
    start<-paste0(year.i, "-01-01")
    end<-paste0(year.i+1,"-01-01")
    imagecoll<-ee$ImageCollection('LANDSAT/LC08/C01/T1_8DAY_EVI')$filterDate(start,end)
    image  <- imagecoll$select('EVI')$toBands()

    #12. Extract buffer mean EVI values----
    image.evi <- image$reduceRegions(collection=data.buff,
                                         reducer=ee$Reducer$mean(),
                                         scale=30)

    #13. Export EVI task to google drive----
    task_vector <- ee_table_to_drive(collection=image.evi,
                                    description=paste0("EVI_hr_",year.i,"-",j),
                                    folder="MCP3",
                                    timePref=FALSE)
    task_vector$start()
#    ee_monitoring(task_vector) # optional

    #14. Load ALAN image collection----
    start<-paste0(year.i, "-01-01")
    end<-paste0(year.i,"-12-31")
    alancoll<-ee$ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG')$filterDate(start,end)$select("avg_rad")$toBands()

    #15. Extract buffer mean ALAN values----
    image.alan <- alancoll$reduceRegions(collection=data.buff,
                                     reducer=ee$Reducer$mean(),
                                     scale=500)

    #16. Export ALAN task to google drive----
    task_vector <- ee_table_to_drive(collection=image.alan,
                                     description=paste0("ALAN_hr_",year.i,"-",j),
                                     folder="MCP3",
                                     timePref=FALSE)
    task_vector$start()
    #    ee_monitoring(task_vector) # optional

    #17. Load human modification image----
    hmi <- ee$ImageCollection('CSP/HM/GlobalHumanModification')$toBands()

    #18. Extract human modification point values----
    image.hmi <- hmi$reduceRegions(collection=data.buff,
                                         reducer=ee$Reducer$mean(),
                                         scale=1000)

    #19. Export ALAN task to google drive----
    task_vector <- ee_table_to_drive(collection=image.hmi,
                                     description=paste0("HMI_hr_",year.i,"-",j),
                                     folder="MCP3",
                                     timePref=FALSE)
    task_vector$start()

    #20. Load Dynamic World images----
    if(year.i %in% c(2015:2017)){
      start<-paste0("2015-06-23")
    }
    else{
      start<-paste0(year.i-2,"-01-01")
    }
    end<-paste0(year.i+2,"-12-31")
    lc <- ee$ImageCollection('GOOGLE/DYNAMICWORLD/V1')$filterDate(start,end)

    #21. Select bands of interest and stack----
    crop <- lc$select("crops")$mean()
    water <- lc$select("water")$mean()
    tree <- lc$select("trees")$mean()
    dw <- crop$addBands(water)$addBands(tree)

    #22. Extract DW meanpoint values----
    image.dw <- dw$reduceRegions(collection=data.buff,
                                   reducer=ee$Reducer$mean(),
                                   scale=10)

    #23. Export DW task to google drive----
    task_vector <- ee_table_to_drive(collection=image.dw,
                                     description=paste0("DW_hr_",year.i,"-",j),
                                     folder="MCP4",
                                     timePref=FALSE)
    task_vector$start()
    #ee_monitoring(task_vector, max_attempts=100) # optional

    end_time <- Sys.time()
    
    print(paste0("Finished loop ", j, " of ", loops, " for year ", year.i, " in ", end_time - start_time, " seconds"))
    
  }

}

#24. Get list of GEE output files----

#NEED TO DOWNLOAD GEE RESULTS AND PUT THEM IN WORKING DIRECTORY HERE####
#Couldn't get ee_drive_to_local to work - Error in if (nrow(files_gd) > 0) { : argument is of length zero
files <- data.frame(file = list.files("Interim/GEE/")) |> 
  separate(file, into=c("image", "scale", "year", "loop", "filename"), remove=FALSE) |>
  dplyr::select(file, image, scale) |> 
  mutate(filepath = paste0("Interim/GEE/", file))

files.evi <- files |> 
  dplyr::filter(image=="EVI", 
                scale=="hr") 

files.alan <- files |> 
  dplyr::filter(image=="ALAN", 
                scale=="hr")

files.hmi <- files |> 
  dplyr::filter(image=="HMI", 
                scale=="hr")

files.dw <- files |>
  dplyr::filter(image=="DW",
                scale=="hr")

#25. Read in ones that don't require date matching----
data.hmi <- readr::read_csv(files.hmi$filepath, show_col_types = FALSE) |> 
  dplyr::select(-'system:index', -'.geo', -'date_millis', -Date) |> 
  group_by(ptID) |> 
  mutate(ptIDn = paste0(ptID,"-",row_number())) |> 
  ungroup() |> 
  rename(hmi = mean)

data.dw <- readr::read_csv(files.dw$filepath, show_col_types = FALSE) |>
  dplyr::select(-'system:index', -'.geo', -'date_millis', -Date) |>
  group_by(ptID) |>
  mutate(ptIDn = paste0(ptID,"-",row_number())) |>
  ungroup() |>
  rename(tree = trees, crop = crops)

#26. Import EVI data, pivot, remove nas, and filter to closest temporal match, randomly pick one if two closest dates----
data.evi <- data.frame()
for(i in 1:nrow(files.evi)){
  data.file <- read.csv(files.evi$filepath[i]) |> 
    dplyr::select(-'system.index', -'.geo', -'date_millis', -Date) |> 
    pivot_longer(cols=c(1:46),
                 names_to="eviimagedate",
                 values_to="evi") |> 
    dplyr::filter(!is.na(evi)) |> 
    mutate(eviimagedate = ymd(str_sub(eviimagedate, 2, 9)),
           date = ymd(str_sub(timestamp, 1, 10)),
           evidatediff = as.numeric(abs(date - eviimagedate))) |> 
    group_by(ptIDn) |> 
    mutate(evimindiff = min(evidatediff)) |> 
    dplyr::filter(evidatediff == evimindiff) |>
    sample_n(1) |> 
    ungroup() |> 
    unique()
  
  data.evi <- rbind(data.evi, data.file)
  
}

#27. Import ALAN data, pivot, remove nas, and filter to closest temporal match, randomly pick one if two closest dates----
data.alan <- data.frame()
for(i in 1:nrow(files.alan)){
  data.file <- read.csv(files.alan$filepath[i]) |> 
    dplyr::select(-'system.index', -'.geo', -'date_millis', -Date) |> 
    pivot_longer(cols=c(1:12),
                 names_to="alanimagedate",
                 values_to="alan") |> 
    dplyr::filter(!is.na(alan)) |> 
    mutate(alanimagedate = ymd(str_sub(alanimagedate, 2, 9)),
           date = ymd(str_sub(timestamp, 1, 10)),
           alandatediff = as.numeric(abs(date - alanimagedate))) |> 
    group_by(ptIDn) |> 
    mutate(alanmindiff = min(alandatediff, na.rm=TRUE)) |> 
    dplyr::filter(alandatediff == alanmindiff) |>
    sample_n(1) |> 
    ungroup() |> 
    unique()
  
  data.alan <- rbind(data.alan, data.file)
  
}

#28. Put all data sources together----
data.join <- data.evi |> 
  mutate(timestamp = ymd_hms(timestamp)) |> 
  full_join(data.hmi) |> 
  full_join(data.dw) |> 
  full_join(data.alan |> 
              mutate(timestamp = ymd_hms(timestamp)) |> 
              dplyr::select(-date)) |> 
  full_join(trackingdata |> 
               mutate(timestamp = ymd_hms(timestamp)) |> 
              dplyr::select(ptIDn, X, Y)) |> 
  unique()

#29. Read in tracking data for season info----
load("Interim/CONIMCP_Habitat.Rdata")
dat.hab <- dat.country |> 
  dplyr::filter(Type != "Band") |> 
  group_by(PinpointID) |> 
  mutate(row=row_number()) |> 
  ungroup() |> 
  mutate(ptID = paste0(PinpointID,"-", row)) |> 
  arrange(ptID) |> 
  dplyr::select(ptID, Sex, Year, Season, Winter)

#30. Put everything together, filter out winter migration points----
data.covs <- data.join |> 
  left_join(dat.hab) |> 
  separate(ptID, into=c("PinpointID", "n"), remove=FALSE) |> 
  dplyr::select(PinpointID, Sex, ptID, Radius, Type, timestamp, Season, Winter, X, Y, evidatediff, alandatediff, evi, alan, hmi, tree, crop, water) |> 
  dplyr::filter(!is.na(tree),
                !is.na(evi),
                !is.na(alan),
                !is.na(hmi),
                !Season=="WinterMig")

#31. Take out points with less than n available points with covs or without a used point----
n <- 20

pt.n.0 <- table(data.covs$ptID, data.covs$Type) |> 
  data.frame() |> 
  rename(ptID=Var1, Type=Var2) |> 
  dplyr::filter(Type=="Available",
                Freq < n) 
table(pt.n.0$Freq)
nrow(pt.n.0)

pt.n.1 <- table(data.covs$ptID, data.covs$Type) |> 
  data.frame() |> 
  rename(ptID=Var1, Type=Var2) |> 
  dplyr::filter(Type=="Used",
                Freq < 1) 
table(pt.n.1$Freq)

data.n <- data.covs |> 
  dplyr::filter(!ptID %in% pt.n.0$ptID,
                !ptID %in% pt.n.1$ptID)

#32. Randomly sample to 25 available points per point----
set.seed(1234)
data.sub <- data.n |> 
  dplyr::filter(Type=="Available") |> 
  group_by(ptID) |> 
  sample_n(n) |> 
  ungroup() |> 
  rbind(data.n |> 
              dplyr::filter(Type=="Used")) |> 
  mutate(used = ifelse(Type=="Used", 1, 0))

#33. Write to csv----
write.csv(data.sub, "Interim/Covariates_landscape.csv", row.names=FALSE)
