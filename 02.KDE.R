library(tidyverse)
library(lubridate)
library(sf)
library(sp)
library(adehabitatHR)
library(data.table)
library(MuMIn)
library(gridExtra)

options(scipen = 999)

#1. Read in data----
ids.exclude <- read.csv("Data/ExclusionYearForBirdsWith2BreedingYears.csv")

dat.hab <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  mutate(PtYear = year(ymd(Date))) %>% 
  anti_join(ids.exclude)

dat.mig <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/DataPaper/coni_migration_data/01_PinPoint2217_StopoverClassified.csv") %>% 
  dplyr::filter(stopover==1,
                count >= 5) %>% 
  mutate(ID = paste0("2217-", season, "-", cluster))

#2. Identify individuals that don't have enough points for KDE ----
dat.use <- dat.hab   %>% 
  dplyr::filter(Season %in% c("Breed", "Winter")) %>% 
  dplyr::filter(Type != "Band", 
                sun==1)

ids <- data.frame(table(dat.use$PinpointID, dat.use$Season, dat.use$Winter)) %>% 
  dplyr::filter(Freq>5) %>% 
  rename(PinpointID=Var1, Season=Var2, Winter=Var3) %>% 
  mutate(PinpointID = as.integer(as.character(PinpointID)),
         Winter = as.integer(as.character(Winter)))

#3. Separate data into KDE and not and add migration stopover points----
#remove banding points for these
dat.kde <- dat.use %>% 
  inner_join(ids) %>% 
  mutate(ID=paste0(PinpointID, "-", Season, "-", Winter)) %>% 
  dplyr::select(ID, PinpointID, Population, Mass, Wing, Sex, Winter, Long, Lat)

#4. Reproject to UTM
dat.kde.m <- dat.kde %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates()

#5. Calculate KDE for individuals with >= 5 points (minimum required)----
dat.kde.sp <- SpatialPointsDataFrame(coords=dat.kde.m, 
                                   data=data.frame(dat.kde$ID),
                                   proj4string = CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "))

kd <- kernelUD(dat.kde.sp, grid = 1000, extent=2, h="href", same4all=FALSE)

#6. Calculate area of 95% isopleth----
kd.area<- kernel.area(kd, percent=95,
                           unin="m", unout="km2",
                           standardize = FALSE) %>% 
  data.frame() %>% 
  transpose(keep.names="ID") %>% 
  mutate(ID=str_sub(ID, 2, 30)) %>% 
  rename(HRarea=V1) %>% 
    separate(ID, into=c("PinpointID", "Season", "Winter")) %>% 
  mutate(PinpointID = as.integer(PinpointID),
         Winter=as.integer(Winter))

ggplot(kd.area) +
  geom_histogram(aes(x=HRarea)) +
  facet_wrap(~Season, scales="free")

#7. Join to data----
dat.area <- kd.area %>%
  right_join(dat.kde) %>% 
  group_by(PinpointID) %>% 
  mutate(count=n()) %>% 
  ungroup() %>% 
  dplyr::select(PinpointID, Population, Season, Winter, count, Sex, Mass, Wing, HRarea) %>% 
  unique() %>% 
  mutate(radius = sqrt(HRarea/3.1416)) %>% 
  left_join(dat.hab %>% 
              dplyr::filter(Season %in% c("Breed", "Winter")) %>% 
              group_by(PinpointID, Winter) %>% 
              summarize(Lat=mean(Lat),
                        Long=mean(Long)) %>% 
              ungroup())

table(dat.area$PinpointID)

#write.csv(dat.area, "KDEArea.csv", row.names = FALSE)

#8. Effects on home range size----
#Visualize
ggplot(dat.area) +
  geom_point(aes(x=count, y=HRarea, colour=factor(Sex))) +
  facet_wrap(~Season, scales="free") #Seems like maybe a relationship for breeding but not for winter

ggplot(dat.area %>% dplyr::filter(Wing > 0)) +
  geom_point(aes(x=Wing, y=HRarea, colour=factor(Sex))) +
  facet_wrap(~Season, scales="free")

ggplot(dat.area) +
  geom_point(aes(x=Lat, y=HRarea, colour=factor(Sex))) +
  facet_wrap(~Season, scales="free")

ggplot(dat.area) +
  geom_point(aes(x=Long, y=HRarea, colour=factor(Sex))) +
  facet_wrap(~Season, scales="free")

#Test for effects of season, sex, sample size
lm.season <- lm(HRarea ~ Season*count + Season*Sex, data=dat.area, na.action="na.fail")
dredge(lm.season)
summary(lm.season)

lm.area.use <- lm(HRarea ~ Sex*season, data=dat.area.lm)
summary(lm.area.use)
newdat <- data.frame(expand.grid(Sex=unique(dat.area.lm$Sex), season=unique(dat.area.lm$season)))
pred.area <- data.frame(predict(lm.area.use, newdat, se=TRUE)) %>% 
  cbind(newdat)

ggplot(pred.area) +
#  geom_ribbon(aes(x=count, ymin=fit-1.96*se.fit, ymax=fit+1.96*se.fit, group=season), fill="grey85") +
  geom_point(aes(x=Sex, y=fit, colour=season))
#  geom_point(aes(x=Sex, y=HRarea, colour=season), data=dat.area.lm)

#10. Calculate KDE for stopover data----
dat.mig.m <- dat.mig %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates()

dat.mig.sp <- SpatialPointsDataFrame(coords=dat.mig.m, 
                                     data=data.frame(dat.mig$ID),
                                     proj4string = CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "))

kd <- kernelUD(dat.mig.sp, grid = 1000, extent=2, h="href", same4all=FALSE)

#11. Calculate area of 95% isopleth for stopover data----
kd.mig.area<- kernel.area(kd, percent=95,
                      unin="m", unout="km2",
                      standardize = FALSE) %>% 
  data.frame() %>% 
  transpose(keep.names="ID") %>% 
  mutate(ID=str_sub(ID, 2, 30)) %>% 
  rename(HRarea=V1) %>% 
  separate(ID, into=c("PinpointID", "season", "cluster")) %>% 
  mutate(PinpointID = as.integer(PinpointID),
         cluster=as.integer(cluster))

#12. Put together and calculate means----
dat.area.all <- kd.mig.area %>% 
  left_join(dat.mig) %>% 
  group_by(PinpointID, season, cluster) %>% 
  summarize(Lat = mean(Latitude),
            Long = mean(Longitude),
            HRarea = mean(HRarea)) %>% 
  ungroup() %>% 
  rename(Season = season,
         id = cluster) %>% 
  mutate(Sex = "M") %>% 
  rbind(dat.area %>% 
          rename(id = Winter) %>% 
          dplyr::select(PinpointID, Season, id, Lat, Long, HRarea, Sex))  %>% 
  mutate(migseason = Season,
         Season = ifelse(Season %in% c("fallmig", "springmig"), "migration", Season))


area <- dat.area.all %>% 
  group_by(Season, Sex) %>% 
  summarize(area.mean = mean(HRarea),
            area.sd = sd(HRarea)) %>% 
  ungroup() %>% 
  rbind(dat.area.all %>% 
          group_by(Season) %>% 
          summarize(area.mean = mean(HRarea),
                    area.sd = sd(HRarea)) %>% 
          ungroup() %>% 
          mutate(Sex = NA)) %>% 
  mutate(radius.mean = sqrt(area.mean/3.1416),
         radius.sd = sqrt(area.sd/3.1416))
area

write.csv(area, "KDEAreaMean.csv", row.names = FALSE)

#TO DO: NEED TO PULL OUT SAMPLE SIZES FOR THESE####

#13. Save out one individual for an example as shapefiles----
dat.kde.i <- dat.kde %>% 
  dplyr::filter(PinpointID==825)

dat.kde.m <- dat.kde.i %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates()

dat.kde.sp <- SpatialPointsDataFrame(coords=dat.kde.m, 
                                     data=data.frame(dat.kde.i$ID),
                                     proj4string = CRS("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "))

kd <- kernelUD(dat.kde.sp, grid = 1000, extent=2, h="href", same4all=FALSE)

kd.shp.95 <- getverticeshr(kd, 95) %>% 
  st_as_sf() %>% 
  mutate(iso=95)
kd.shp.75 <- getverticeshr(kd, 75) %>% 
  st_as_sf() %>% 
  mutate(iso=75)
kd.shp.50 <- getverticeshr(kd, 50) %>% 
  st_as_sf() %>% 
  mutate(iso=50)
kd.shp.25 <- getverticeshr(kd, 25) %>% 
  st_as_sf() %>% 
  mutate(iso=25)
kd.shp.05 <- getverticeshr(kd, 5) %>% 
  st_as_sf() %>% 
  mutate(iso=5)

kd.shp <- rbind(kd.shp.95, kd.shp.75, kd.shp.50, kd.shp.25, kd.shp.05) %>% 
  separate(id, into=c("PinpointID", "Season"), remove=FALSE)

kd.shp.breed <- kd.shp %>% 
  dplyr::filter(Season=="Breed2")
dat.kde.i.breed <- dat.kde.i %>% 
  dplyr::filter(Season=="Breed2") %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(dat.kde.i %>% 
          dplyr::filter(Season=="Breed2"))

kd.shp.winter <- kd.shp %>% 
  dplyr::filter(Season=="Winter2")
dat.kde.i.winter <- dat.kde.i %>% 
  dplyr::filter(Season=="Winter2") %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(dat.kde.i %>% 
          dplyr::filter(Season=="Winter2"))

breed <- ggplot() +
  geom_sf(data=kd.shp.breed, aes(fill=iso)) +
  geom_point(data=dat.kde.i.breed, aes(x=X, y=Y, colour=doy)) +
  scale_colour_viridis_c()

winter <- ggplot() +
  geom_sf(data=kd.shp.winter, aes(fill=iso)) +
  geom_point(data=dat.kde.i.winter, aes(x=X, y=Y, colour=doy)) +
  scale_colour_viridis_c()

grid.arrange(breed, winter)

write_sf(kd.shp, "Shapefiles/ExampleKDE.shp")
write.csv(dat.kde.i, "Shapefiles/ExampleKDEData.csv", row.names = FALSE)
