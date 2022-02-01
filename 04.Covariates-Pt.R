#https://www.mdpi.com/2072-4292/13/20/4154

library(sf)
library(tidyverse)
library(lubridate)
library(rgee)
library(data.table)
library(usdm)

#1. Initialize rgee----
#ee_install()
ee_Initialize()
ee_check()

#2. Settings----

#Number of samples
n <- 25

#Distance kernel
kern <- ee$Kernel$euclidean(radius=10000, units="meters")

#2. Write functions----
#2a. Function to add property with time in milliseconds
add_date<-function(feature) {
  date <- ee$Date(ee$String(feature$get("Date")))$millis()
  feature$set(list(date_millis=date))
}

#3. Import data----
trackingdata <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting_pt.csv", header = T) %>% 
  rename(timestamp=DateTime) %>% 
  arrange(ptID) %>% 
  fill(timestamp, .direction="down") %>% 
  group_by(ptID) %>% 
  mutate(ptIDn = paste0(ptID,"-",row_number())) %>% 
  ungroup()
head(trackingdata)

#4. Wrangle data----
trackingdata$Date <- as.POSIXct(trackingdata$timestamp, format = "%Y-%m-%d %H:%M:%S", tz="UTC") 
trackingdata$Date <- as.factor(trackingdata$Date)
trackingdata$year <- year(ymd_hms(trackingdata$timestamp))
trackingdata$Date <- sub(" ", "T", trackingdata$Date) #Put in a format that can be read by javascript

#Sample a few points for testing----
#trackingdata <- trackingdata %>% 
#  sample_n(1000)

#5. Set up loop to go through each year----
years <- unique(trackingdata$year)

data.year <- list()
data.out <- list()

for(i in 1:length(years)){
  
  start_time <- Sys.time()
  
  year.i <- years[i]
  
  data.i <- trackingdata %>% 
    dplyr::filter(year==year.i) %>% 
    arrange(X, Y) %>% 
    mutate(row = row_number(),
           n = ceiling(row/1000))
  
  loops <- max(data.i$n)
  
  for(j in 1:loops){
    
    data.j <- data.i %>% 
      dplyr::filter(n==j)
    
    #8. Create sf object----
    datasf <- st_as_sf(data.j, coords = c('X','Y'), crs = 3857) %>% 
      st_transform(crs=4326)
    
    #9. Send data to GEE----
    data <- sf_as_ee(datasf)
    
    #10. Transform day into milliseconds----
    data <-data$map(add_date)
    
    #11. Load EVI image collection----
    start<-paste0(year.i, "-01-01")
    end<-paste0(year.i,"-12-31")
    imagecoll<-ee$ImageCollection('LANDSAT/LC08/C01/T1_8DAY_EVI')$filterDate(start,end)
    
    #12. Extract EVI point values----
    data.evi <- ee_extract(
      x = imagecoll,
      y = datasf,
      scale = 30,
      sf = FALSE
    )
    
    #13. Pivot, remove nas, and filter to closest temporal match, randomly pick one if two closest dates----
    set.seed(1234)
    data.piv <- data.evi %>% 
      pivot_longer(cols=c(10:ncol(data.evi)),
                   names_to="imagedate",
                   values_to="EVI") %>% 
      dplyr::filter(!is.na(EVI)) %>% 
      mutate(imagedate = ymd(str_sub(imagedate, 2, 9)),
             timestamp = ymd(str_sub(timestamp, 1, 10)),
             datediff = as.numeric(abs(timestamp - imagedate))) %>% 
      group_by(ptIDn) %>% 
      mutate(mindiff = min(datediff)) %>% 
      dplyr::filter(datediff == mindiff) %>%
      sample_n(1) %>% 
      ungroup() %>% 
      mutate(year = year.i) %>% 
      dplyr::select(-timestamp) %>% 
      right_join(data.j)
    
    #14. Load Copernicus image----
    lc <- ee$Image('COPERNICUS/Landcover/100m/Proba-V-C3/Global/2017')
    
    #15. Extract Copernicus point values----
    data.lc <- ee_extract(
      x = lc,
      y = datasf,
      scale = 100,
      sf = FALSE
    )
    
    #16. Select only the discrete classification band----
    lcdc <- lc$select('discrete_classification')
    
    #17. Create raster of distance to open water----
    water <- lcdc$mask(lcdc$eq(80))
    openwaterdist <- water$distance(kern, skipMasked=FALSE)$rename("openwaterdist")
    
    #18. Get point value of distance to open water----
    data.openwater <- ee_extract(
      x = openwaterdist,
      y = datasf,
      scale = 100,
      sf = FALSE
    )
    
    #19. Create raster of distance to wetland----
    wetland <- lcdc$mask(lcdc$eq(90))
    wetlanddist <- wetland$distance(kern, skipMasked=FALSE)$rename("wetlanddist")
    
    #20. Get point value of distance to open water----
    data.wetland <- ee_extract(
      x = wetlanddist,
      y = datasf,
      scale = 100,
      sf = FALSE
    )
    
    #21. Create raster of distance to crop----
    crop <- lcdc$mask(lcdc$eq(40))
    cropdist <- crop$distance(kern, skipMasked=FALSE)$rename("cropdist")
    
    #22. Get point value of distance to open water----
    data.crop <- ee_extract(
      x = cropdist,
      y = datasf,
      scale = 100,
      sf = FALSE
    )
  
    #23. Put all data sources together----
    data.cov <- full_join(data.piv, data.lc) %>% 
      full_join(data.openwater) %>% 
      full_join(data.wetland) %>% 
      full_join(data.crop)

    #24. Save to list----
    data.year[[j]] <- data.cov
    
    end_time <- Sys.time()
    
    print(paste0("Finished loop ", j, " of ", loops, " for year ", year.i, " in ", end_time - start_time, " minutes"))
    
  }
  
  data.out[[i]] <- rbindlist(data.year)
  
}

#18. Convert to dataframe & save----
data.all <- rbindlist(data.out, fill=TRUE) %>% 
  mutate(doy = yday(ymd_hms(timestamp))) %>% 
  right_join(trackingdata) %>% 
  unique()

write.csv(data.all, "Data/Covariates_pt_raw.csv")

data.n <- table(data.all$ptIDn) %>% 
  data.frame() %>% 
  rename(ptIDn = Var1) %>% 
  dplyr::filter(Freq > 1) %>% 
  left_join(data.all)

#19. Read in tracking data for season info----
dat.hab <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  dplyr::filter(Type != "Band") %>% 
  group_by(PinpointID) %>% 
  mutate(row=row_number()) %>% 
  ungroup() %>% 
  mutate(ptID = paste0(PinpointID,"-", row)) %>% 
  arrange(ptID) %>% 
  dplyr::select(ptID, Year, Season, Winter)

#20. Put everything together, filter out winter migration points----
data.covs <- data.all %>% 
  rename_with(~gsub(pattern=".coverfraction", replacement="", .x)) %>% 
  rename_with(~gsub(pattern="water.", replacement="", .x)) %>% 
  mutate(water = permanent + seasonal) %>% 
  left_join(dat.hab) %>% 
  separate(ptID, into=c("pinpointID", "n"), remove=FALSE) %>% 
dplyr::select(pinpointID, ptID, Radius, Type, timestamp, Season, Winter, X, Y, datediff, EVI, bare, crops, grass, moss, shrub, tree, water) %>% 
  dplyr::filter(!is.na(tree),
                !is.na(EVI),
                !Season=="WinterMig")

#21. Take out points with less than 20 available points with covs----
pt.n.0 <- table(data.covs$ptID, data.covs$Type) %>% 
  data.frame() %>% 
  rename(ptID=Var1, Type=Var2) %>% 
  dplyr::filter(Type=="Available",
                Freq < n) 
table(pt.n.0$Freq)

pt.n.1 <- table(data.covs$ptID, data.covs$Type) %>% 
  data.frame() %>% 
  rename(ptID=Var1, Type=Var2) %>% 
  dplyr::filter(Type=="Used",
                Freq < 1) 

data.n <- data.covs %>% 
  dplyr::filter(!ptID %in% pt.n.0$ptID,
                !ptID %in% pt.n.1$ptID)

#22. Randomly sample to 25 available points per point----
set.seed(1234)
data.sub <- data.n %>% 
  dplyr::filter(Type=="Used") %>% 
  rbind(data.n %>% 
          dplyr::filter(Type=="Available") %>% 
          group_by(ptID) %>% 
          sample_n(n) %>% 
          ungroup()) %>% 
  mutate(used = ifelse(Type=="Used", 1, 0))

#23. Split into seasons, add BirdID field----
season <- unique(data.sub$Season)

data.season <- data.frame()
for(i in 1:length(season)){
  
  data.season.i <- data.sub %>% 
    dplyr::filter(Season==season[i]) %>% 
    dplyr::select(pinpointID) %>% 
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
           evi.s = scale(EVI))
  
  data.season <- rbind(data.season, data.season.i)
  
}

#24. Write to csv----
write.csv(data.season, "Data/Covariates_pt.csv", row.names=FALSE)

#25. VIF----
data.season <- read.csv("Data/Covariates_pt.csv")

data.vif <- data.season %>% 
  dplyr::select(tree.s, grass.s, shrub.s, bare.s, crops.s, water.s, evi.s)
cor(data.vif)
#grass and tree
vif(data.vif)

data.vif <- data.season %>% 
  dplyr::select(tree.s, shrub.s, bare.s, crops.s, water.s, evi.s)
cor(data.vif)
#grass and tree
vif(data.vif)
