library(sf)
library(tidyverse)
library(lubridate)
library(rgee)
library(data.table)

#1. Initialize rgee----
ee_Initialize()
ee_check()

#2. Settings---
#Buffer radius
rad <- 200

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
data.raw <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting_5km.csv", header = T)

trackingdata <- data.raw %>% 
  st_as_sf(coords = c('X','Y'), crs = 3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  cbind(data.raw %>% 
          dplyr::select(-X, -Y)) %>% 
  rename(timestamp=DateTime) %>% 
  arrange(ptID) %>% 
  fill(timestamp, .direction="down") %>% 
  group_by(ptID) %>% 
  mutate(ptIDn = paste0(ptID, "-", row_number())) %>% 
  ungroup()
  
#5. Wrangle data----
trackingdata$Date <- as.POSIXct(trackingdata$timestamp, format = "%Y-%m-%d %H:%M:%S", tz="UTC") 
trackingdata$Date <- as.factor(trackingdata$Date)
trackingdata$year <- year(ymd_hms(trackingdata$timestamp))
trackingdata$Date <- sub(" ", "T", trackingdata$Date) #Put in a format that can be read by javascript

#Sample a few points for testing----
#trackingdata <- trackingdata %>% sample_n(1000)

#6. Set up loop to go through each year----
years <- unique(trackingdata$year)

for(i in 1:length(years)){
  
  start_time <- Sys.time()
  
  year.i <- years[i]
  
  data.i <- trackingdata %>% 
    dplyr::filter(year==year.i) %>% 
    mutate(row = row_number(),
           n = ceiling(row/1000))
  
  loops <- max(data.i$n)
  
  for(j in 1:loops){
    
    data.j <- data.i %>% 
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
                                    folder="MCP",
                                    timePref=FALSE)
    task_vector$start()
#    ee_monitoring(task_vector) # optional
    
    #14. Load Copernicus image----
    lc <- ee$Image('COPERNICUS/Landcover/100m/Proba-V-C3/Global/2017')
    
    #15. Extract buffer mean lc values----
    image.lc <- lc$reduceRegions(collection=data.buff, 
                                    reducer=ee$Reducer$mean(), 
                                    scale=30)
    
    #16. Export Copernicus task to google drive----
    task_vector <- ee_table_to_drive(collection=image.lc,
                                     description=paste0("LC_hr_",year.i,"-",j),
                                     folder="MCP",
                                     timePref=FALSE)
    task_vector$start()
#    ee_monitoring(task_vector) # optional
#    img <- ee_drive_to_local(task = task_vector)
    
    
    end_time <- Sys.time()
    
    print(paste0("Finished loop ", j, " of ", loops, " for year ", year.i, " in ", end_time - start_time, " seconds"))
    
  }
  
}

#NEED TO DOWNLOAD GEE RESULTS AND PUT THEM IN WORKING DIRECTORY HERE####
#Couldn't get ee_drive_to_local to work - Error in if (nrow(files_gd) > 0) { : argument is of length zero

#17. Get list of GEE output files----
files <- data.frame(file = list.files("Data/GEE/")) %>% 
  separate(file, into=c("image", "scale", "year", "loop", "filename"), remove=FALSE) %>%
  dplyr::select(file, image, scale) %>% 
  mutate(filepath = paste0("Data/GEE/", file))

files.evi <- files %>% 
  dplyr::filter(image=="EVI", 
                scale=="hr") 

files.lc <- files %>% 
  dplyr::filter(image=="LC", 
                scale=="hr")

#18. Import Copernicus data----
data.lc <- readr::read_csv(files.lc$filepath, show_col_types = FALSE) %>% 
  dplyr::select(-'system:index', -'.geo', -'date_millis', -Date) %>% 
  group_by(ptID) %>% 
  mutate(ptIDn = paste0(ptID,"-",row_number())) %>% 
  ungroup()

#19. Import EVI data, pivot, remove nas, and filter to closest temporal match, randomly pick one if two closest dates----
data.evi <- data.frame()
for(i in 1:nrow(files.evi)){
  data.file <- read.csv(files.evi$filepath[i]) %>% 
    dplyr::select(-'system.index', -'.geo', -'date_millis', -Date) %>% 
    pivot_longer(cols=c(1:46),
                 names_to="imagedate",
                 values_to="evi") %>% 
    dplyr::filter(!is.na(evi)) %>% 
    mutate(imagedate = ymd(str_sub(imagedate, 2, 9)),
           date = ymd(str_sub(timestamp, 1, 10)),
           datediff = as.numeric(abs(date - imagedate))) %>% 
    group_by(ptIDn) %>% 
    mutate(mindiff = min(datediff)) %>% 
    dplyr::filter(datediff == mindiff) %>%
    sample_n(1) %>% 
    ungroup() %>% 
    mutate(year = year.i) %>% 
    unique()
  
  data.evi <- rbind(data.evi, data.file)
  
}

#20. Put two data sources together----
data.all <- data.evi %>% 
  mutate(timestamp = ymd_hms(timestamp)) %>% 
  full_join(data.lc) %>% 
  full_join(trackingdata %>% 
               mutate(timestamp = ymd_hms(timestamp)) %>% 
              dplyr::select(ptIDn, X, Y)) %>% 
  unique()

#21. Read in tracking data for season info----
dat.hab <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  dplyr::filter(Type != "Band") %>% 
  group_by(PinpointID) %>% 
  mutate(row=row_number()) %>% 
  ungroup() %>% 
  mutate(ptID = paste0(PinpointID,"-", row)) %>% 
  arrange(ptID) %>% 
  dplyr::select(PinpointID, ptID, Year, Season, Winter)

#22. Put everything together, filter out winter migration points----
data.covs <- data.all %>% 
  rename_with(~gsub(pattern=".coverfraction", replacement="", .x)) %>% 
  rename_with(~gsub(pattern="water-", replacement="", .x)) %>% 
  mutate(water = permanent + seasonal) %>% 
  left_join(dat.hab) %>% 
  dplyr::select(PinpointID, ptID, Radius, Type, timestamp, Season, Winter, datediff, evi, bare, crops, grass, moss, shrub, tree, water) %>% 
  dplyr::filter(!is.na(tree),
                !is.na(evi),
                !Season=="WinterMig")

#23. Take out IDs with less than 20 available points with covs & IDs with no used point----
pt.n.0 <- table(data.covs$ptID, data.covs$Type) %>% 
  data.frame() %>% 
  rename(ptID=Var1, Type=Var2) %>% 
  dplyr::filter(Type=="Available",
                Freq < 20) 
table(pt.n.0$Freq)

pt.n.1 <- table(data.covs$ptID, data.covs$Type) %>% 
  data.frame() %>% 
  rename(ptID=Var1, Type=Var2) %>% 
  dplyr::filter(Type=="Used",
                Freq < 1) 

data.n <- data.covs %>% 
  dplyr::filter(!ptID %in% pt.n.0$ptID,
                !ptID %in% pt.n.1$ptID)

#24. Randomly sample to 20 available points per point----
set.seed(1234)
data.sub <- data.n %>% 
  dplyr::filter(Type=="Available") %>% 
  group_by(ptID) %>% 
  sample_n(20) %>% 
  ungroup() %>% 
  rbind(data.n %>% 
              dplyr::filter(Type=="Used")) %>% 
  mutate(used = ifelse(Type=="Used", 1, 0))

#25. Split into seasons, add BirdID field----
season <- unique(data.sub$Season)

data.season <- data.frame()
for(i in 1:length(season)){
  
  data.season.i <- data.sub %>% 
    dplyr::filter(Season==season[i]) %>% 
    dplyr::select(PinpointID) %>% 
    unique() %>% 
    mutate(BirdID = row_number()) %>% 
    right_join(data.sub %>% 
                 dplyr::filter(Season==season[i])) %>% 
    arrange(BirdID, ptID, used) %>% 
    mutate(tree.s = scale(tree),
           grass.s = scale(grass),
           shrub.s = scale(shrub),
           bare.s = scale(bare),
           crops.s = scale(crops),
           water.s = scale(water),
           evi.s = scale(evi))
  
  data.season <- rbind(data.season, data.season.i)
  
}

#24. Write to csv----
write.csv(data.season, "Data/Covariates_hr.csv", row.names=FALSE)