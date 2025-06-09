# ---
# title: CONI roost habitat selection - selecting available points for habitat selection
# author: Elly Knight
# created: November 9, 2021
# updated: June 9, 2025
# ---

# Preamble - load packages
library(tidyverse) # data wrangling
library(sf) # working with shps
library(data.table) # list handling

#1. Read in used data----
load("Interim/CONIMCP_habitat.Rdata")

#wrangle - filter to roosting only, no banding locations, no winter migration locations, no nesting females
dat.hab <- dat.country |> 
  dplyr::filter(Type != "Band",
                Season !="WinterMig",
                behaviour=="roost",
                !(Sex=="F" & Season=="Breed")) |> 
  group_by(PinpointID) |> 
  mutate(row=row_number()) |> 
  ungroup() |> 
  mutate(ptID = paste0(PinpointID,"-", row)) |> 
  arrange(ptID)

#2. Point level: Random points within 300 m----
n.pt <- 30

#make the buffer
buff.pt <- dat.hab |> 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) |> 
  st_transform(crs=3857) |> 
  st_buffer(dist=300) |> 
  arrange(PinpointID, Season)

IDs.pt <- data.frame(ptID = rep(buff.pt$ptID, n.pt)) |> 
  arrange(ptID)

#pick the available points
set.seed(1234)
pts.pt <- st_sample(buff.pt, size=rep(n.pt, nrow(buff.pt))) |> 
  st_coordinates() |> 
  data.frame() |> 
  cbind(IDs.pt) |> 
  mutate(DateTime = NA,
         Type="Available",
         Radius="pt")

#put it together
dat.pt <- dat.hab |> 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) |> 
  st_transform(crs=3857) |> 
  st_coordinates() |> 
  cbind(dat.hab) |> 
  mutate(Type="Used",
         Radius="pt") |> 
  dplyr::select(ptID, X, Y, DateTime, Type, Radius) |> 
  rbind(pts.pt) 

#make it a shp
dat.pt.sf <- dat.pt |> 
  st_as_sf(coords=c("X", "Y"), crs=3857) |> 
  cbind(dat.pt |> 
          dplyr::select(X, Y))

#save
write.csv(dat.pt, "Interim/CONIMCP_Local.csv", row.names = FALSE)
write_sf(dat.pt.sf, "Shapefiles/CONIMCP_Local.shp", append=FALSE)

#4. Home range level: Random points within HR size----
#decide the radius
area <- read.csv("Interim/KDEAreaMean.csv") |> 
  rename(Season = Season2) |> 
  mutate(radius = round(radius.mean)) |> 
  dplyr::select(Season, radius)
area #use 1, 5, 50

n.hr <- 30

#some wrangling
buff.hr.breed <- dat.hab |> 
  dplyr::filter(Season=="Breed") |> 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) |> 
  st_transform(crs=3857) |> 
  st_buffer(dist=5000) |> 
  arrange(PinpointID)

buff.hr.winter <- dat.hab |> 
  dplyr::filter(Season=="Winter") |> 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) |> 
  st_transform(crs=3857) |> 
  st_buffer(dist=1000) |> 
  arrange(PinpointID)

buff.hr.migration <- dat.hab |> 
  dplyr::filter(Season %in% c("FallMig", "SpringMig")) |> 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) |> 
  st_transform(crs=3857) |> 
  st_buffer(dist=50000) |> 
  arrange(PinpointID)

buff.hr <- rbind(buff.hr.breed, buff.hr.winter, buff.hr.migration)

IDs.hr <- data.frame(ptID = rep(buff.hr$ptID, n.hr)) |> 
  arrange(ptID)

#pick the available points
set.seed(1234)
pts.hr <- st_sample(buff.hr, size=rep(n.hr, nrow(buff.hr))) |> 
  st_coordinates() |> 
  data.frame() |> 
  cbind(IDs.hr) |> 
  mutate(DateTime = NA,
         Type="Available",
         Radius="hr")

#put back together
dat.hr <- dat.hab |> 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) |> 
  st_transform(crs=3857) |> 
  st_coordinates() |> 
  cbind(dat.hab) |> 
  mutate(Type="Used",
         Radius="hr") |> 
  dplyr::select(ptID, X, Y, DateTime, Type, Radius) |> 
  rbind(pts.hr)

#make a shp
dat.hr.sf <- dat.hr |> 
  st_as_sf(coords=c("X", "Y"), crs=3857) |> 
  cbind(dat.hr |> 
          dplyr::select(X, Y))

#save
write.csv(dat.hr, "Interim/CONIMCP_Landscape.csv", row.names = FALSE)
write_sf(dat.hr.sf, "Shapefiles/CONIMCP_Landscape.shp", append=FALSE)