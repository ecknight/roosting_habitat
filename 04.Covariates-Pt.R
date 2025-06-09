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
kern <- ee$Kernel$euclidean(radius=20000, units="meters")

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
trackingdata$month <- month(ymd_hms(trackingdata$timestamp))
trackingdata$Date <- sub(" ", "T", trackingdata$Date) #Put in a format that can be read by javascript

#Sample a few points for testing----
#trackingdata <- trackingdata %>% 
#  sample_n(1000)

#5. Set up loop to go through each year----
years <- unique(sort(trackingdata$year))

data.out <- list()

for(i in 1:length(years)){
  
  year.i <- years[i]
  
  data.i <- trackingdata %>% 
    dplyr::filter(year==year.i) %>% 
    arrange(X, Y) %>% 
    mutate(row = row_number(),
           n = ceiling(row/1000))
  
  loops <- max(data.i$n)
  data.year <- list()
  for(j in 1:loops){
    
    start_time <- Sys.time()
    
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
    evicoll<-ee$ImageCollection('LANDSAT/LC08/C01/T1_8DAY_EVI')$filterDate(start,end)
    
    #12. Extract EVI point values----
    data.evi <- ee_extract(
      x = evicoll,
      y = datasf,
      scale = 30,
      sf = FALSE
    )
    
    #13. Pivot, remove nas, and filter to closest temporal match, randomly pick one if two closest dates----
    set.seed(1234)
    data.evi.p <- data.evi %>% 
      pivot_longer(cols=c(10:ncol(data.evi)),
                   names_to="eviimagedate",
                   values_to="EVI") %>% 
      dplyr::filter(!is.na(EVI)) %>% 
      mutate(eviimagedate = ymd(str_sub(eviimagedate, 2, 9)),
             timestamp = ymd(str_sub(timestamp, 1, 10)),
             evidatediff = as.numeric(abs(timestamp - eviimagedate))) %>% 
      group_by(ptIDn) %>% 
      mutate(evimindiff = min(evidatediff, na.rm=TRUE)) %>% 
      dplyr::filter(evidatediff == evimindiff) %>%
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
    
    # #17. Create raster of distance to open water----
    # water <- lcdc$mask(lcdc$eq(80))
    # openwaterdist <- water$distance(kern, skipMasked=FALSE)$rename("openwaterdist")
    # 
    # #18. Get point value of distance to open water----
    # data.openwater <- ee_extract(
    #   x = openwaterdist,
    #   y = datasf,
    #   scale = 100,
    #   sf = FALSE
    # )
    # 
    # #19. Create raster of distance to wetland----
    # wetland <- lcdc$mask(lcdc$eq(90))
    # wetlanddist <- wetland$distance(kern, skipMasked=FALSE)$rename("wetlanddist")
    # 
    # #20. Get point value of distance to wetland----
    # data.wetland <- ee_extract(
    #   x = wetlanddist,
    #   y = datasf,
    #   scale = 100,
    #   sf = FALSE
    # )
    # 
    # #21. Create raster of distance to crop----
    # crop <- lcdc$mask(lcdc$eq(40))
    # cropdist <- crop$distance(kern, skipMasked=FALSE)$rename("cropdist")
    # 
    # #22. Get point value of distance to crop----
    # data.crop <- ee_extract(
    #   x = cropdist,
    #   y = datasf,
    #   scale = 100,
    #   sf = FALSE
    # )
    
    #28. Load ALAN image collection----
    start<-paste0(year.i, "-01-01")
    end<-paste0(year.i,"-12-31")
    alancoll<-ee$ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG')$filterDate(start,end)$select("avg_rad")
    
    #29. Extract ALAN point values----
    data.alan <- ee_extract(
      x = alancoll,
      y = datasf,
      scale = 500,
      sf = FALSE
    )
    
    #30. Pivot, remove nas, and filter to closest temporal match, randomly pick one if two closest dates----
    set.seed(1234)
    data.alan.p <- data.alan %>% 
      pivot_longer(cols=c(10:ncol(data.alan)),
                   names_to="alanimagedate",
                   values_to="ALAN") %>% 
      dplyr::filter(!is.na(ALAN)) %>% 
      mutate(alanimagedate = ymd(str_sub(alanimagedate, 2, 9)),
             timestamp = ymd(str_sub(timestamp, 1, 10)),
             alandatediff = as.numeric(abs(timestamp - alanimagedate))) %>% 
      group_by(ptIDn) %>% 
      mutate(alanmindiff = min(alandatediff, na.rm=TRUE)) %>% 
      dplyr::filter(alandatediff == alanmindiff) %>%
      sample_n(1) %>% 
      ungroup() %>% 
      mutate(year = year.i) %>% 
      dplyr::select(-timestamp) %>% 
      right_join(data.j)
    
    #31. Load human modification image----
    hmi <- ee$ImageCollection('CSP/HM/GlobalHumanModification')
    
    #32. Extract human modification point values----
    data.hmi <- ee_extract(
      x = hmi,
      y = datasf,
      scale = 1000,
      sf = FALSE
    )
    
    #33. Load Dynamic World images----
    if(year.i %in% c(2015:2017)){
      start<-paste0("2015-06-23")
    }
    else {start<-paste0(year.i-2,"-01-01")}
    end<-paste0(year.i+2,"-12-31")
    lc <- ee$ImageCollection('GOOGLE/DYNAMICWORLD/V1')$filterDate(start,end)
    
    #34. Select bands of interest----
    crop <- lc$select("crops")$mean()
    water <- lc$select("water")$mean()
    tree <- lc$select("trees")$mean()
    
    #35. Extract dynamic world point values----
    data.tree.dw <- ee_extract(
      x = tree,
      y = datasf,
      scale = 10,
      sf = FALSE
    )
    
    data.crop.dw <- ee_extract(
      x = crop,
      y = datasf,
      scale = 10,
      sf = FALSE
    )
    
    data.water.dw <- ee_extract(
      x = water,
      y = datasf,
      scale = 10,
      sf = FALSE
    )
    
    # #25. Mask out forest cover less than 30%----
    # cover <- tree$gte(0.3)$selfMask()
    # 
    # #26. Calculate patch size----
    # patch <- cover$connectedComponents(connectedness=ee$Kernel$plus(1), maxSize=500)
    # patchid <- patch$select('labels')
    # patchn <- patchid$connectedPixelCount(maxSize=500, eightConnected=FALSE)
    # nsize <- tree$pixelArea()
    # patchsize <- patchn$multiply(nsize)
    # 
    # #27. Get point value of patch size----
    # data.patch <- try(ee_extract(
    #   x = patchsize,
    #   y = datasf,
    #   scale = 10,
    #   sf = FALSE
    # ))
    # 
    # if(ncol(data.patch)==11){
    #   data.patch <- data.patch %>% 
    #     rename(patch = labels) %>% 
    #     mutate(patch = ifelse(is.na(patch), 0, patch))
    # }
    # else{
    #   data.patch <- data.patch %>% 
    #     mutate(patch=0)
    # }

    #36. Put all data sources together----

    
    data.cov <- full_join(data.evi.p, data.lc) %>% 
      # full_join(data.openwater) %>% 
      # full_join(data.wetland) %>% 
      # full_join(data.crop) %>% 
      # full_join(data.patch) %>% 
      full_join(data.alan.p) %>% 
      full_join(data.hmi) %>% 
      full_join(data.tree.dw) %>%
      full_join(data.crop.dw) %>%
      full_join(data.water.dw)

    #37. Save to list----
    data.year[[j]] <- data.cov
    
    end_time <- Sys.time()
    
    print(paste0("Finished loop ", j, " of ", loops, " for year ", year.i, " in ", end_time - start_time, " minutes"))
    
  }
  
  data.out[[i]] <- rbindlist(data.year, fill=TRUE)
  
}

#38. Convert to dataframe & save----
data.all <- rbindlist(data.out, fill=TRUE) %>% 
  mutate(doy = yday(ymd_hms(timestamp))) %>% 
  right_join(trackingdata) %>% 
  unique() %>% 
  rename(HMI = X2016_gHM,
         treedw = trees,
         cropdw = crops, 
         waterdw = water)

write.csv(data.all, "Data/Covariates_pt_raw.csv")

data.all <- read.csv("Data/Covariates_pt_raw.csv")

data.n <- table(data.all$ptIDn) %>% 
  data.frame() %>% 
  rename(ptIDn = Var1) %>% 
  dplyr::filter(Freq > 1) %>% 
  left_join(data.all)

#39. Read in tracking data for season info----
dat.hab <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  dplyr::filter(Type != "Band") %>% 
  group_by(PinpointID) %>% 
  mutate(row=row_number()) %>% 
  ungroup() %>% 
  mutate(ptID = paste0(PinpointID,"-", row)) %>% 
  arrange(ptID) %>% 
  dplyr::select(ptID, Sex, Year, Season, Winter)

#40. Put everything together, filter out winter migration points----
data.covs <- data.all %>% 
  rename_with(~gsub(pattern=".coverfraction", replacement="", .x)) %>% 
  mutate(water = water.permanent + water.seasonal,
         # openwaterdist = ifelse(is.na(openwaterdist), 30000, openwaterdist),
         # wetlanddist = ifelse(is.na(wetlanddist), 30000, wetlanddist),
         # cropdist = ifelse(is.na(cropdist), 30000, cropdist),
         #waterdist1 = ifelse(openwaterdist<wetlanddist, 1, 0),
         # waterdist = ifelse(waterdist1==1, openwaterdist, wetlanddist)
         )%>% 
  left_join(dat.hab) %>% 
  separate(ptID, into=c("pinpointID", "n"), remove=FALSE) %>% 
dplyr::select(pinpointID, Sex, ptID, Radius, Type, timestamp, Season, Winter, X, Y, evidatediff, alandatediff, EVI, bare, crops, grass, moss, shrub, tree, water,
#              waterdist, cropdist,
              patch, ALAN, HMI, treedw, cropdw, waterdw) %>% 
  dplyr::filter(!is.na(tree),
                !is.na(EVI),
                !is.na(treedw),
                !Season=="WinterMig")

#41. Take out points with less than n available points with covs----
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

#42. Randomly sample to n available points per point----
set.seed(1234)
data.sub <- data.n %>% 
  dplyr::filter(Type=="Used") %>% 
  rbind(data.n %>% 
          dplyr::filter(Type=="Available") %>% 
          group_by(ptID) %>% 
          sample_n(n) %>% 
          ungroup()) %>% 
  mutate(used = ifelse(Type=="Used", 1, 0))

#44. Write to csv----
write.csv(data.sub, "Data/Covariates_pt.csv", row.names=FALSE)

#45. VIF----
data.out <- read.csv("Data/Covariates_pt.csv")

data.vif <- data.out %>% 
  dplyr::select(EVI, ALAN, HMI, treedw, cropdw, waterdw)
cor(data.vif)
vif(data.vif)
