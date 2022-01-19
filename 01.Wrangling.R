library(tidyverse)
library(suncalc)
library(lubridate)
library(data.table)
library(adehabitatLT)
library(meanShiftR) #for density cluster analysis
library(sf)
library(rworldmap)
library(readxl)

options(scipen = 999)

whemi <- map_data("world", region=c("Canada", 
                                    "USA", 
                                    "Mexico",
                                    "Guatemala", 
                                    "Belize", 
                                    "El Salvador",
                                    "Honduras", 
                                    "Nicaragua", 
                                    "Costa Rica",
                                    "Panama", 
                                    "Jamaica", 
                                    "Cuba", 
                                    "The Bahamas",
                                    "Haiti", 
                                    "Dominican Republic", 
                                    "Antigua and Barbuda",
                                    "Dominica", 
                                    "Saint Lucia", 
                                    "Saint Vincent and the Grenadines", 
                                    "Barbados",
                                    "Grenada",
                                    "Trinidad and Tobago",
                                    "Colombia",
                                    "Venezuela",
                                    "Guyana",
                                    "Suriname",
                                    "Ecuador",
                                    "Peru",
                                    "Brazil",
                                    "Bolivia",
                                    "Paraguay",
                                    "Chile",
                                    "Argentina",
                                    "Uruguay")) %>% 
  dplyr::filter(!group%in%c(258:264))

#1. Read in data----
dat <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv") %>% 
  mutate(Date = as.Date(str_sub(DateTime, 1, 10)),
         DateTime = ymd_hms(paste0(Date, Time))) %>% 
  dplyr::select(PinpointID, Population, Mass, Wing, Sex, Type, DateTime, Date, Time, Year, Lat, Long, BandDist, WintDist, GCD, Season, Season2, Winter) %>% 
  unique()
  
#2. Calculate sun times----
dat.sun <- getSunlightTimes(data=dat%>% 
                               rename(lat = Lat, lon = Long, date=Date) %>% 
                               dplyr::select(date, lat, lon)) %>% 
  dplyr::select(sunrise, sunset) %>% 
  cbind(dat) %>% 
  mutate(sun = ifelse(DateTime > sunrise & DateTime < sunset, 1, 0))

#3. Calculate NSD----
traj <- as.ltraj(xy=dat.sun[,c("Long", "Lat")],
                 id=dat.sun$PinpointID,
                 date=dat.sun$DateTime,
                 proj4string = CRS("+proj=longlat +datum=WGS84")) 

dat.traj <- rbindlist(traj) %>% 
  dplyr::select(dist, R2n, abs.angle, rel.angle) %>% 
  cbind(dat.sun) %>% 
  mutate(doy=yday(DateTime)) %>% 
  group_by(PinpointID) %>% 
  mutate(R2n2 = max(R2n) - R2n) %>% 
  ungroup()

#4. Prepare for clustering----
dat.coord <- dat.traj %>% 
  st_as_sf(crs=4326, coords=c("Long", "Lat")) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(dat.traj)

ids <- dat.coord %>% 
  dplyr::select(PinpointID) %>% 
  unique()

#5. Mean shift classification----
dat.shift <- data.frame()
for(i in 1:nrow(ids)){
  
  dat.i <- dat.coord %>% 
    dplyr::filter(PinpointID==ids$PinpointID[i])
  
  mat1 <- matrix(dat.i$X)
  mat2 <- matrix(dat.i$Y)
  mat <- cbind(mat1, mat2)
  
  shift <- meanShift(mat,
                         algorithm="KDTREE",
                         bandwidth=c(1,1))
  
  dat.shift <- dat.i %>% 
    mutate(cluster = shift[[1]]) %>% 
    rbind(dat.shift)
}

dat.clust <- dat.shift %>% 
  group_by(PinpointID, cluster) %>% 
  mutate(count=n()) %>% 
  ungroup() %>% 
  mutate(stopover = ifelse(count > 3, 1, 0))

#6. Clean wintering data----
dat.wint <- dat.clust %>% 
  dplyr::filter(stopover==1) %>% 
  dplyr::filter(Lat < 25) %>% 
  group_by(PinpointID, cluster) %>% 
  mutate(minDate = min(DateTime),
         maxDate = max(DateTime),
         duration = maxDate - minDate) %>% 
  ungroup() %>% 
  dplyr::filter(duration > 21) %>% 
  mutate(Season="Winter") %>% 
  dplyr::select(-Winter)

#Visualize
plot.wint <- ggplot() +
  geom_point(data=dat.wint, aes(x=Long, y=Lat, colour=doy), size=3, alpha=0.7) +
  labs(x = "", y = "") +
  theme_bw() +
  scale_colour_viridis_c() +
  facet_wrap(~PinpointID, ncol=10, scales="free")

ggsave(plot.wint, file="Figures/MeanShift_wint.jpeg", width=20, height=10, unit="in")

#7. individuals with 2 wintering grounds
ids.wint <- dat.wint %>% 
  dplyr::select(PinpointID, Season, cluster) %>% 
  unique() %>% 
  group_by(PinpointID, Season) %>% 
  summarize(n=n()) %>% 
  dplyr::filter(n==2) %>% 
  left_join(dat.wint %>% 
              dplyr::select(PinpointID, Season, cluster, minDate) %>% 
              unique()) %>% 
  mutate(Winter = ifelse(minDate==min(minDate), 1, 2)) %>% 
  ungroup()

dat.wint2 <- dat.wint %>% 
  left_join(ids.wint) %>% 
  mutate(Winter = ifelse(is.na(Winter), 1, Winter)) %>% 
  dplyr::filter(!(PinpointID==489 & Lat > -8.5 & Type!="Band"),
                !(PinpointID==480 & Long < -67 & Type!="Band"),
                !(PinpointID==443 & Long < -28.89)) %>% 
  dplyr::select(PinpointID, Population, Mass, Wing, Sex, Type, DateTime, Date, doy, Time, sun, Year, Lat, Long, BandDist, WintDist, GCD, Season, Winter, dist, R2n, abs.angle, rel.angle, cluster, count)

#8. Winter migration----
dat.wintmig <-  dat.clust %>% 
  dplyr::filter(Lat < 25) %>% 
  group_by(PinpointID, cluster) %>% 
  mutate(minDate = min(DateTime),
         maxDate = max(DateTime),
         duration = maxDate - minDate) %>% 
  ungroup() %>% 
  dplyr::filter(duration > 21) %>% 
  dplyr::filter(PinpointID %in% ids.wint$PinpointID,
                stopover==0) %>% 
  mutate(Season="WinterMig",
         Winter=0) %>% 
  dplyr::select(PinpointID, Population, Mass, Wing, Sex, Type, DateTime, Date, doy, Time, sun, Year, Lat, Long, BandDist, WintDist, GCD, Season, Winter, dist, R2n, abs.angle, rel.angle, cluster, count)

#9. Clean breeding data----
dat.ids <- dat.wint %>% 
  dplyr::select(PinpointID) %>% 
  unique()

dat.breed <- dat.traj %>% 
  dplyr::filter(PinpointID %in% dat.ids$PinpointID) %>% 
  left_join(dat.clust) %>% 
  dplyr::filter(Lat > 25,
                stopover==1 | Type=="Band") %>% 
  group_by(PinpointID) %>% 
  mutate(maxpt = ifelse(Lat==max(Lat), cluster, 0),
         maxclust = max(maxpt),
         useclust = ifelse(cluster==maxclust, 1, 0)) %>% 
  mutate(TypeNo = ifelse(Type=="Band", 0, 1),
         TypeSum = sum(TypeNo),
         BandUse = case_when(TypeSum==0 ~ 1,
                             is.na(maxclust) & Type=="Band" ~ 1,
                             PinpointID %in% c(480, 489) ~1)) %>% 
  ungroup() %>% 
  mutate(Season="Breed",
         Winter=0)  %>% 
  dplyr::filter((useclust==1 & Type != "Band") | BandUse==1) %>% 
  dplyr::filter(!(PinpointID==443 & Lat < 46)) %>% 
  dplyr::select(PinpointID, Population, Mass, Wing, Sex, Type, DateTime, Date, doy, Time, sun, Year, Lat, Long, BandDist, WintDist, GCD, Season, Winter, dist, R2n, abs.angle, rel.angle, cluster, count)

#Visualize
plot.breed <- ggplot() +
  geom_point(data=dat.breed, aes(x=Long, y=Lat, colour=factor(Type)), size=3, alpha=0.7) +
  labs(x = "", y = "") +
  theme_bw() +
  scale_colour_viridis_d() +
  facet_wrap(~PinpointID, ncol=10, scales="free")

ggsave(plot.breed, file="Figures/MeanShift_Breed.jpeg", width=20, height=10, unit="in")

#10. Fall migration----
dat.fall <- dat.clust %>% 
  dplyr::select(PinpointID, Population, Mass, Wing, Sex, Type, DateTime, Date, doy, Time, sun, Year, Lat, Long, BandDist, WintDist, GCD,  dist, R2n, abs.angle, rel.angle, cluster, count, Season2) %>% 
  anti_join(dat.wint) %>% 
  anti_join(dat.wintmig) %>% 
  anti_join(dat.breed) %>% 
  dplyr::filter(year(DateTime)==Year,
                !Type=="Band",
                !Season2=="WinterMig") %>% 
  mutate(Season="FallMig",
         Winter=0) %>% 
  dplyr::select(PinpointID, Population, Mass, Wing, Sex, Type, DateTime, Date, doy, Time, sun, Year, Lat, Long, BandDist, WintDist, GCD, Season, Winter, dist, R2n, abs.angle, rel.angle, cluster, count)

#11. Spring migration----
dat.spring <- dat.clust %>% 
  dplyr::select(PinpointID, Population, Mass, Wing, Sex, Type, DateTime, Date, doy, Time, sun, Year, Lat, Long, BandDist, WintDist, GCD,  dist, R2n, abs.angle, rel.angle, cluster, count, Season2) %>% 
  anti_join(dat.wint) %>% 
  anti_join(dat.wintmig) %>% 
  anti_join(dat.breed) %>% 
  dplyr::filter(year(DateTime)>Year,
                !Type=="Band",
                !Season2=="WinterMig") %>% 
  mutate(Season="SpringMig",
         Winter=0) %>% 
  dplyr::select(PinpointID, Population, Mass, Wing, Sex, Type, DateTime, Date, doy, Time, sun, Year, Lat, Long, BandDist, WintDist, GCD, Season, Winter, dist, R2n, abs.angle, rel.angle, cluster, count)

#12. All data----
dat.all <- rbind(dat.wint2, dat.wintmig, dat.breed, dat.fall, dat.spring) %>% 
  dplyr::select(PinpointID, Population, Mass, Wing, Sex, Type, DateTime, Year, sun, Lat, Long, Season, Winter, cluster, count) %>% 
  mutate(cluster = as.numeric(cluster))

table(dat.all$Season, dat.all$PinpointID)

#13. Add archival pinpoint tag data----
band <- read_excel("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/tbl_band.xlsx") %>% 
  dplyr::filter(PinpointID=="2217")

dat.2217 <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/DataPaper/coni_migration_data/01_PinPoint2217_StopoverClassified.csv") %>% 
  rename(Lat = Latitude, Long = Longitude, Season = season, DateTime = date) %>% 
  mutate(PinpointID=2217,
         Population = 8,
         Type="Retrieval",
         DateTime = ymd_hms(DateTime),
         sunrise= ymd_hms(sunrise),
         sunset = ymd_hms(sunset),
         Year = year(DateTime),
         sun = ifelse(DateTime > sunrise & DateTime < sunset, 1, 0),
         Winter = ifelse(Season=="Winter", 1, 0),
         Mass = band$Mass, 
         Wing = band$Wing,
         Sex = band$Sex,
         Season = case_when(Season=="breed" ~ "Breed",
                            Season=="winter" ~ "Winter",
                            Season=="fallmig" ~ "FallMig",
                            Season=="springmig" ~ "SpringMig")) %>% 
    dplyr::select(PinpointID, Population, Mass, Wing, Sex, Type, DateTime, Year, sun, Lat, Long, Season, Winter, cluster, count)

dat.all.2217 <- rbind(dat.all, dat.2217)

write.csv(dat.all.2217, "Data/CONIMCP_CleanDataAll_Habitat.csv", row.names = FALSE)

#14. Just roosting data---
dat.roost <- dat.all.2217 %>% 
  dplyr::filter(sun==1 | Type=="Band")

table(dat.roost$Season, dat.roost$PinpointID)

#Visualize
plot.shift <- ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_point(data=dat.roost, aes(x=Long, y=Lat, colour=factor(Season)), size=3, alpha=0.7) +
  labs(x = "", y = "") +
  xlim(c(-170, -30)) +
  theme_bw() +
  scale_colour_viridis_d() +
  facet_wrap(~PinpointID, ncol=10)

#ggsave(plot.shift, file="Figures/MeanShift_Map.jpeg", width=20, height=10, unit="in")

#15. Assign country to each point & remove points over the gulf----
countries <- getMap(resolution='low')
points <- dat.roost %>% 
  dplyr::select(Long, Lat) %>% 
  SpatialPoints(proj4string = CRS(proj4string(countries)))
points.country <- over(points, countries)

dat.country <- dat.roost %>% 
  mutate(country = as.character(points.country$ADMIN),
         region = as.character(points.country$continent)) %>% 
  dplyr::filter(!is.na(country)) %>% 
  arrange(PinpointID, DateTime)

write.csv(dat.country, "Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv", row.names = FALSE)

#16. Final numebrs----
length(unique(dat.country$PinpointID))
length(unique(dat.country$Population))
nrow(dat.country)
