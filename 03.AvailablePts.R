library(tidyverse)
library(sf)
library(data.table)

#1. Read in used data----
dat.hab <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  dplyr::filter(Type != "Band") %>% 
  group_by(PinpointID) %>% 
  mutate(row=row_number()) %>% 
  ungroup() %>% 
  mutate(ID = paste0(PinpointID,"-", row)) %>% 
  arrange(ID)

#3. Point level: Random points within 200 m----
n.pt <- 100

IDs <- data.frame(ID = rep(dat.hab$ID, n.pt)) %>% 
  arrange(ID)

buff.pt <- dat.hab %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_buffer(dist=200) %>% 
  arrange(PinpointID, Season)

set.seed(1234)
pts.pt <- st_sample(buff.pt, size=rep(n.pt, nrow(buff.pt))) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(IDs) %>% 
  mutate(DateTime = NA,
         Type="Available",
         Radius="200m")

dat.pt <- dat.hab %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(dat.hab) %>% 
  mutate(Type="Used",
         Radius="200m") %>% 
  dplyr::select(X, Y, ID, DateTime, Type, Radius) %>% 
  rbind(pts.pt)
  
write.csv(dat.pt, "Data/CONIMCP_CleanDataAll_Habitat_Roosting_200m.csv")

dat.pt.sf <- dat.pt %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  cbind(dat.pt %>% 
          dplyr::select(X, Y))

plot(dat.pt.sf)

write_sf(dat.pt.sf, "Shapefiles/CONIMCP_CleanDataAll_Habitat_Roosting_200m.shp", append=FALSE)

#4. Home range level: Random points within 5 km----
n.hr <- 20

IDs <- data.frame(ID = rep(dat.hab$ID, n.hr)) %>% 
  arrange(ID)

buff.pt <- dat.hab %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_buffer(dist=5000) %>% 
  arrange(PinpointID, Season)

set.seed(1234)
pts.pt <- st_sample(buff.pt, size=rep(n.hr, nrow(buff.pt))) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(IDs) %>% 
  mutate(DateTime = NA,
         Type="Available",
         Radius="5km")

dat.pt <- dat.hab %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(dat.hab) %>% 
  mutate(Type="Used",
         Radius="5km") %>% 
  dplyr::select(X, Y, ID, DateTime, Type, Radius) %>% 
  rbind(pts.pt)

write.csv(dat.pt, "Data/CONIMCP_CleanDataAll_Habitat_Roosting_5km.csv")

dat.pt.sf <- dat.pt %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  cbind(dat.pt %>% 
          dplyr::select(X, Y))

plot(dat.pt.sf)

write_sf(dat.pt.sf, "Shapefiles/CONIMCP_CleanDataAll_Habitat_Roosting_5km.shp", append=FALSE)

#5. Landscape level: Random points within 100 km----
n.land <- 20

IDs <- data.frame(ID = rep(dat.hab$ID, n.land)) %>% 
  arrange(ID)

buff.pt <- dat.hab %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_buffer(dist=100000) %>% 
  arrange(PinpointID, Season)

set.seed(1234)
pts.pt <- st_sample(buff.pt, size=rep(n.land, nrow(buff.pt))) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(IDs) %>% 
  mutate(DateTime = NA,
         Type="Available",
         Radius="100km")

dat.pt <- dat.hab %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(dat.hab) %>% 
  mutate(Type="Used",
         Radius="100km") %>% 
  dplyr::select(X, Y, ID, DateTime, Type, Radius) %>% 
  rbind(pts.pt)

write.csv(dat.pt, "Data/CONIMCP_CleanDataAll_Habitat_Roosting_100km.csv")

dat.pt.sf <- dat.pt %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  cbind(dat.pt %>% 
          dplyr::select(X, Y))

#plot(dat.pt.sf)

write_sf(dat.pt.sf, "Shapefiles/CONIMCP_CleanDataAll_Habitat_Roosting_100km.shp", append=FALSE)
