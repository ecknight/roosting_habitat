library(tidyverse)
library(sf)
library(jagsUI)
library(coda)
library(data.table)

options(scipen = 9999)

#1. Read in tracking data & put together----
raw <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  dplyr::filter(Type != "Band") %>% 
  group_by(PinpointID) %>% 
  mutate(row=row_number()) %>% 
  ungroup() %>% 
  mutate(ptID = paste0(PinpointID,"-", row)) %>% 
  dplyr::select(ptID, country, region, GCD, Sex, Mass, Population, Lat, Long, Season)

#3. Calculate distance to gulf----
gom <- read_sf("Shapefiles/GulfOfMexico.shp") %>% 
  st_make_valid()

dat.sf <- raw %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326)

dat.dist <- data.frame(distance = as.numeric(st_distance(dat.sf, gom))/1000) %>% 
  cbind(raw) %>% 
  mutate(position = if_else(Lat > 25, "north", "south")) %>% 
  dplyr::select(ptID, country, region, GCD, Sex, Mass, Population, Season, distance, position, Lat, Long) %>% 
  unique() %>% 
  mutate(distance.rd = round(distance, -3),
         distance.rd = ifelse(distance.rd > 4000, 4000, distance.rd))

#Visualize
ggplot() +
  geom_sf(data=gom) +
  geom_point(data=dat.dist, aes(x=Long, y=Lat, colour=distance)) +
  facet_grid(Season ~ position) +
  scale_colour_viridis_c()

#4. Put together with covariate data----
pt <- read.csv("Data/Covariates_pt.csv") %>% 
  rename(PinpointID = pinpointID, 
         evi = EVI)

hr <- read.csv("Data/Covariates_hr.csv") %>% 
  st_as_sf(coords=c("X", "Y"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(read.csv("Data/Covariates_hr.csv") %>% 
          dplyr::select(-X, -Y))

land <- read.csv("Data/Covariates_land.csv") %>% 
  st_as_sf(coords=c("X", "Y"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(read.csv("Data/Covariates_land.csv") %>% 
          dplyr::select(-X, -Y))

dat <- rbind(pt, hr, land) %>% 
  mutate(scale = case_when(Radius=="200m" ~ "pt",
                           Radius=="5km" ~ "hr",
                           Radius=="100km" ~ "land"),
         scale = factor(scale, levels=c("pt", "hr", "land"))) %>% 
  mutate(water2 = water^2,
         water2.s = water.s^2) %>% 
  left_join(dat.dist) %>% 
  mutate(distance.rd.dir = case_when(Season=="FallMig" & position=="north" ~ -distance.rd,
                                     Season=="SpringMig" & position=="south" ~ -distance.rd,
                                     !is.na(distance.rd) ~ distance.rd))

#5. Visualize some things----
ggplot(dat %>% dplyr::filter(Season %in% c("FallMig", "SpringMig"))) +
  geom_smooth(aes(x=evi, y=used, colour=Season)) +
  facet_grid(scale ~ distance.rd, scales="free")
