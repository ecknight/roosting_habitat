library(tidyverse)
library(sf)
library(data.table)

#1. Read in used data----
dat.hab <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  dplyr::filter(Type != "Band",
                Season !="WinterMig") %>% 
  group_by(PinpointID) %>% 
  mutate(row=row_number()) %>% 
  ungroup() %>% 
  mutate(ptID = paste0(PinpointID,"-", row)) %>% 
  arrange(ptID)

#2. Point level: Random points within 300 m----
n.pt <- 30

buff.pt <- dat.hab %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_buffer(dist=300) %>% 
  arrange(PinpointID, Season)

IDs.pt <- data.frame(ptID = rep(buff.pt$ptID, n.pt)) %>% 
  arrange(ptID)

set.seed(1234)
pts.pt <- st_sample(buff.pt, size=rep(n.pt, nrow(buff.pt))) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(IDs.pt) %>% 
  mutate(DateTime = NA,
         Type="Available",
         Radius="pt")

dat.pt <- dat.hab %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(dat.hab) %>% 
  mutate(Type="Used",
         Radius="pt") %>% 
  dplyr::select(ptID, X, Y, DateTime, Type, Radius) %>% 
  rbind(pts.pt) 
  
write.csv(dat.pt, "Data/CONIMCP_CleanDataAll_Habitat_Roosting_pt.csv", row.names = FALSE)

dat.pt.sf <- dat.pt %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  cbind(dat.pt %>% 
          dplyr::select(X, Y))

#plot(dat.pt.sf)

write_sf(dat.pt.sf, "Shapefiles/CONIMCP_CleanDataAll_Habitat_Roosting_pt.shp", append=FALSE)

#4. Home range level: Random points within HR size----
area <- read.csv("KDEAreaMean.csv") %>% 
  dplyr::filter(is.na(Sex)) %>% 
  mutate(radius = round(radius.mean)) %>% 
  dplyr::select(Season, radius)
area

n.hr <- 30

buff.hr.breed <- dat.hab %>% 
  dplyr::filter(Season=="Breed") %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_buffer(dist=5000) %>% 
  arrange(PinpointID)

buff.hr.winter <- dat.hab %>% 
  dplyr::filter(Season=="Winter") %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_buffer(dist=1000) %>% 
  arrange(PinpointID)

buff.hr.migration <- dat.hab %>% 
  dplyr::filter(Season %in% c("FallMig", "SpringMig")) %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_buffer(dist=50000) %>% 
  arrange(PinpointID)

buff.hr <- rbind(buff.hr.breed, buff.hr.winter, buff.hr.migration)

IDs.hr <- data.frame(ptID = rep(buff.hr$ptID, n.hr)) %>% 
  arrange(ptID)

set.seed(1234)
pts.hr <- st_sample(buff.hr, size=rep(n.hr, nrow(buff.hr))) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(IDs.hr) %>% 
  mutate(DateTime = NA,
         Type="Available",
         Radius="hr")

dat.hr <- dat.hab %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(dat.hab) %>% 
  mutate(Type="Used",
         Radius="hr") %>% 
  dplyr::select(ptID, X, Y, DateTime, Type, Radius) %>% 
  rbind(pts.hr)

write.csv(dat.hr, "Data/CONIMCP_CleanDataAll_Habitat_Roosting_hr.csv", row.names=FALSE)

dat.hr.sf <- dat.hr %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  cbind(dat.hr %>% 
          dplyr::select(X, Y))

#plot(dat.hr.sf)

write_sf(dat.hr.sf, "Shapefiles/CONIMCP_CleanDataAll_Habitat_Roosting_hr.shp", append=FALSE)