# ---
# title: CONI roost habitat selection - extracting local scale covariate with GEE
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

#2. Write functions----
#2a. Function to add property with time in milliseconds
add_date<-function(feature) {
  date <- ee$Date(ee$String(feature$get("Date")))$millis()
  feature$set(list(date_millis=date))
}

#3. Import data----
trackingdata <- read.csv("Interim/CONIMCP_Local.csv", header = T) |> 
  rename(timestamp=DateTime) |> 
  arrange(ptID) |> 
  fill(timestamp, .direction="down") |> 
  group_by(ptID) |> 
  mutate(ptIDn = paste0(ptID,"-",row_number())) |> 
  ungroup()

#4. Wrangle data----
trackingdata$Date <- as.POSIXct(trackingdata$timestamp, format = "%Y-%m-%d %H:%M:%S", tz="UTC") 
trackingdata$Date <- as.factor(trackingdata$Date)
trackingdata$year <- year(ymd_hms(trackingdata$timestamp))
trackingdata$month <- month(ymd_hms(trackingdata$timestamp))
trackingdata$Date <- sub(" ", "T", trackingdata$Date) #Put in a format that can be read by javascript

#Sample a few points for testing----
#trackingdata <- trackingdata |> 
#  sample_n(1000)

#5. Set up loop to go through each year----
years <- unique(sort(trackingdata$year))

data.out <- list()
for(i in 1:length(years)){
  
  year.i <- years[i]
  
  data.i <- trackingdata |> 
    dplyr::filter(year==year.i) |> 
    arrange(X, Y) |> 
    mutate(row = row_number(),
           n = ceiling(row/1000))
  
  loops <- max(data.i$n)
  data.year <- list()
  for(j in 1:loops){
    
    start_time <- Sys.time()
    
    data.j <- data.i |> 
      dplyr::filter(n==j)
    
    #8. Create sf object----
    datasf <- st_as_sf(data.j, coords = c('X','Y'), crs = 3857) |> 
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
    data.evi.p <- data.evi |> 
      pivot_longer(cols=c(10:ncol(data.evi)),
                   names_to="eviimagedate",
                   values_to="EVI") |> 
      dplyr::filter(!is.na(EVI)) |> 
      mutate(eviimagedate = ymd(str_sub(eviimagedate, 2, 9)),
             timestamp = ymd(str_sub(timestamp, 1, 10)),
             evidatediff = as.numeric(abs(timestamp - eviimagedate))) |> 
      group_by(ptIDn) |> 
      mutate(evimindiff = min(evidatediff, na.rm=TRUE)) |> 
      dplyr::filter(evidatediff == evimindiff) |>
      sample_n(1) |> 
      ungroup() |> 
      mutate(year = year.i) |> 
      dplyr::select(-timestamp) |> 
      right_join(data.j)
    
    #14. Load ALAN image collection----
    start<-paste0(year.i, "-01-01")
    end<-paste0(year.i,"-12-31")
    alancoll<-ee$ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG')$filterDate(start,end)$select("avg_rad")
    
    #15. Extract ALAN point values----
    data.alan <- ee_extract(
      x = alancoll,
      y = datasf,
      scale = 500,
      sf = FALSE
    )
    
    #16. Pivot, remove nas, and filter to closest temporal match, randomly pick one if two closest dates----
    set.seed(1234)
    data.alan.p <- data.alan |> 
      pivot_longer(cols=c(10:ncol(data.alan)),
                   names_to="alanimagedate",
                   values_to="ALAN") |> 
      dplyr::filter(!is.na(ALAN)) |> 
      mutate(alanimagedate = ymd(str_sub(alanimagedate, 2, 9)),
             timestamp = ymd(str_sub(timestamp, 1, 10)),
             alandatediff = as.numeric(abs(timestamp - alanimagedate))) |> 
      group_by(ptIDn) |> 
      mutate(alanmindiff = min(alandatediff, na.rm=TRUE)) |> 
      dplyr::filter(alandatediff == alanmindiff) |>
      sample_n(1) |> 
      ungroup() |> 
      mutate(year = year.i) |> 
      dplyr::select(-timestamp) |> 
      right_join(data.j)
    
    #17. Load human modification image----
    hmi <- ee$ImageCollection('CSP/HM/GlobalHumanModification')
    
    #18. Extract human modification point values----
    data.hmi <- ee_extract(
      x = hmi,
      y = datasf,
      scale = 1000,
      sf = FALSE
    )
    
    #19. Load Dynamic World images----
    if(year.i %in% c(2015:2017)){
      start<-paste0("2015-06-23")
    }
    else {start<-paste0(year.i-2,"-01-01")}
    end<-paste0(year.i+2,"-12-31")
    lc <- ee$ImageCollection('GOOGLE/DYNAMICWORLD/V1')$filterDate(start,end)
    
    #20. Select bands of interest----
    crop <- lc$select("crops")$mean()
    water <- lc$select("water")$mean()
    tree <- lc$select("trees")$mean()
    
    #21. Extract dynamic world point values----
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

    #22. Put all data sources together----
    data.cov <- full_join(data.evi.p, data.alan.p) |> 
      full_join(data.hmi) |> 
      full_join(data.tree.dw) |>
      full_join(data.crop.dw) |>
      full_join(data.water.dw)

    #23. Save to list----
    data.year[[j]] <- data.cov
    
    end_time <- Sys.time()
    
    print(paste0("Finished loop ", j, " of ", loops, " for year ", year.i, " in ", end_time - start_time, " minutes"))
    
  }
  
  data.out[[i]] <- rbindlist(data.year, fill=TRUE)
  
}

#24. Convert to dataframe & save----
data.all <- rbindlist(data.out, fill=TRUE) |> 
  mutate(doy = yday(ymd_hms(timestamp))) |> 
  right_join(trackingdata) |> 
  unique() |> 
  rename(HMI = X2016_gHM,
         tree = trees,
         crop = crops)

write.csv(data.all, "Interim/Covariates_local_raw.csv", row.names = FALSE)

#25. Inventory number of points per used point----
data.all <- read.csv("Interim/Covariates_local_raw.csv")

data.n <- table(data.all$ptID) |> 
  data.frame() |> 
  rename(ptID = Var1) |> 
  dplyr::filter(Freq > 1) |> 
  left_join(data.all)

#26. Read in tracking data for season info----
load("Interim/CONIMCP_Habitat.Rdata")
dat.hab <- dat.country |> 
  dplyr::filter(Type != "Band") |> 
  group_by(PinpointID) |> 
  mutate(row=row_number()) |> 
  ungroup() |> 
  mutate(ptID = paste0(PinpointID,"-", row)) |> 
  arrange(ptID) |> 
  dplyr::select(ptID, Sex, Year, Season, Winter)

#27. Put everything together, filter out winter migration points----
data.covs <- data.all |> 
  left_join(dat.hab) |> 
  separate(ptID, into=c("PinpointID", "n"), remove=FALSE) |> 
  rename(evi = EVI,
         alan = ALAN,
         hmi = HMI) |> 
dplyr::select(PinpointID, Sex, ptID, Radius, Type, timestamp, Season, Winter, X, Y, evidatediff, alandatediff, evi, alan, hmi, tree, crop, water) |> 
  dplyr::filter(!is.na(tree),
                !is.na(evi),
                !is.na(alan),
                !is.na(hmi),
                !Season=="WinterMig")

#28. Take out points with less than n available points with covs or without a used point----
n <- 20

pt.n.0 <- table(data.covs$ptID, data.covs$Type) |> 
  data.frame() |> 
  rename(ptID=Var1, Type=Var2) |> 
  dplyr::filter(Type=="Available",
                Freq < n) 
table(pt.n.0$Freq)

pt.n.1 <- table(data.covs$ptID, data.covs$Type) |> 
  data.frame() |> 
  rename(ptID=Var1, Type=Var2) |> 
  dplyr::filter(Type=="Used",
                Freq < 1) 
table(pt.n.1$Freq)

data.n <- data.covs |> 
  dplyr::filter(!ptID %in% pt.n.0$ptID,
                !ptID %in% pt.n.1$ptID)

#29. Randomly sample to n available points per point----
set.seed(1234)
data.sub <- data.n |> 
  dplyr::filter(Type=="Used") |> 
  rbind(data.n |> 
          dplyr::filter(Type=="Available") |> 
          group_by(ptID) |> 
          sample_n(n) |> 
          ungroup()) |> 
  mutate(used = ifelse(Type=="Used", 1, 0))

#30. Write to csv----
write.csv(data.sub, "Interim/Covariates_local.csv", row.names=FALSE)