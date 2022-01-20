options(scipen = 99999)

library(tidyverse)
library(sf)
library(maps)
library(Cairo)
library(gridExtra)
library(cowplot)
library(raster)
library(lubridate)
library(ggspatial)
library(gganimate)
library(ggmap)
library(ggforce)
library(gridExtra)
library(grid)

my.theme <- theme_classic() +
  theme(text=element_text(size=12, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(size=12, hjust = 0.5))

map.theme <- theme_nothing() +
  theme(text=element_text(size=12, family="Arial"),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.text = element_blank())

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

#https://stackoverflow.com/questions/47749078/how-to-put-a-geom-sf-produced-map-on-top-of-a-ggmap-produced-raster
ggmap_bbox <- function(map) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
  map_bbox <- setNames(unlist(attr(map, "bb")), 
                       c("ymin", "xmin", "ymax", "xmax"))
  
  bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
  
  attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
  attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
  attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
  attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
  map
}

#1. Figure 1 - Study area----
#Deployment metadata
pop <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/tbl_population_abundance.csv") %>% 
  dplyr::filter(Region != "Florida") %>% 
  dplyr::mutate(Region = case_when(Region=="BC coast" ~ "Coastal BC",
                                   Region=="BC Okanagan" ~ "Southcentral BC",
                                   !is.na(Region) ~ as.character(Region))) %>% 
  arrange(-Lat) %>% 
  mutate(order=row_number()) 

#Create breeding ground points with number of individuals deployed
pt.breed <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  mutate(Population = ifelse(Population==6, 7, Population)) %>% 
  dplyr::select(PinpointID, Population) %>% 
  unique() %>%
  left_join(pop) %>% 
  group_by(Population, Lat, Long) %>% 
  summarize(n=n()) %>% 
  mutate(Season="Breeding (# tags)",
         Season = factor(Season))

#All data
dat <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  dplyr::filter(Season!="WinterMig") %>% 
  mutate(DateTime = ymd_hms(DateTime),
         Season = factor(Season, levels=c("Breed", "FallMig", "Winter", "SpringMig"),
                         labels=c("Breeding (# tags)", "Fall Migration", "Winter", "Spring Migration"))) %>% 
  arrange(PinpointID, DateTime) %>% 
  dplyr::filter(Season != "Breeding (# tags)")

studyarea <- ggplot() +
  geom_polygon(data=whemi, aes(x=long, y=lat, group=group), fill="gray70", colour = "gray85", size=0.3) +
  geom_path(data=dat, aes(x=Long, y=Lat, group=PinpointID), colour="grey40", size=0.3) +
  geom_point(data=dat, aes(x=Long, y=Lat, fill=Season), pch=21, colour="grey40") +
  geom_point(data=pt.breed, aes(x=Long, y=Lat), fill="yellow", colour="black", pch=21, size=4) +
  geom_text(data=pt.breed, aes(x=Long, y=Lat, label=n), nudge_y=0, nudge_x=0, size=2.5) +
  xlab("") +
  ylab("") +
  xlim(c(-169, -30)) +
  #  ylim(c(14, 85)) +
  my.theme +
  theme(legend.position = "none",
        plot.margin = unit(c(0,1,-0.5,-0.5), "cm"),
        axis.text.x.bottom = element_text(size=10),
        axis.text.y.left = element_text(size=10)) +
  theme(panel.grid.major = element_line(colour = "gray90")) +
  facet_wrap(~Season)
studyarea

ggsave(studyarea, filename = "Figures/StudyArea.jpeg", height = 6, width = 6)

#2. Figure 2 - KDE & Availability design----
#2a. Histograms----
dat.area.all <- read.csv("KDEArea.csv") %>% 
  mutate(Season = ifelse(Season %in% c("FallMig", "SpringMig"), "Migration", Season),
         Season = ifelse(Season=="Breed", "Breeding", Season),
         Season = factor(Season, levels=c("Breeding", "Winter", "Migration")))

plot.hist <- ggplot(dat.area.all) +
  geom_histogram(aes(x=HRarea)) +
  facet_wrap(~Season, ncol=1, scales="free") +
  xlab("95% isopleth area (km2)") +
  ylab("") +
  my.theme +
  theme(strip.background = element_blank(),
    strip.text.x = element_blank())

#2b. KDE examples----
kd.shp <- read_sf("Shapefiles/ExampleKDE.shp") %>% 
  mutate(iso = factor(iso, levels=c(95, 75, 50, 25, 5))) %>% 
  mutate(Season=ifelse(Season=="SpringMig", "Migration", Season)) %>% 
  st_transform(crs=3857)

dat.kde <- read.csv("Shapefiles/ExampleKDEData.csv") %>% 
  st_as_sf(coords=c("Long", "Lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(read.csv("Shapefiles/ExampleKDEData.csv")) %>% 
  separate(ID, into=c("PinpointID", "Season", "id"), remove=FALSE) %>% 
  mutate(Season=ifelse(Season=="SpringMig", "Migration", Season))

#Breed
#Subset shapefile
breed.shp <- kd.shp %>% 
  dplyr::filter(Season=="Breed")

#Get spatial attributes
center.breed <- dat.kde %>% 
  dplyr::filter(Season=="Breed") %>% 
  summarize(Long = mean(Long),
            Lat = mean(Lat))

center.breed.shp <- breed.shp %>% 
  dplyr::filter(iso==95) %>% 
  st_centroid() %>% 
  st_coordinates() %>% 
  data.frame()

bbox.breed <- st_bbox(breed.shp)

width.breed <- data.frame(width = c(bbox.breed[3] - bbox.breed[1], y = bbox.breed[4] - bbox.breed[2])) %>%
  summarize(max = max(width))

#Get background data
register_google(key="AIzaSyCta9P4x7jGNELznpwlx07VZkkLVk3FP4M")

map.breed <- get_map(center.breed, zoom=11, force=TRUE, maptype="satellite", color="color")
map.breed <- ggmap_bbox(map.breed)

#Plot
plot.kde.breed <- ggmap(map.breed) + 
  coord_sf(crs = st_crs(3857)) +
  geom_sf(data = breed.shp, aes(fill=iso), colour="grey40", alpha = 0.7, inherit.aes=FALSE) +
  geom_point(data=dat.kde %>% dplyr::filter(Season=="Breed"), aes(x=X, y=Y), size=4, pch=21, colour="grey40", fill="grey70") +
  scale_fill_viridis_d(name="Isopleth %", direction=-1) +
  map.theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  xlim(c(center.breed.shp$X-width.breed$max*0.6, center.breed.shp$X+width.breed$max*0.6)) +
  ylim(c(center.breed.shp$Y-width.breed$max*0.6, center.breed.shp$Y+width.breed$max*0.6)) +
  ggsn::scalebar(transform=FALSE,
                 dist=5, dist_unit="km",
                 box.fill=c("grey80", "grey20"),
                 box.color="grey20",
                 location = "bottomleft",
                 height = 0.05,
                 st.dist = 0.025,
                 st.color="grey80",
                 data=breed.shp,
                 anchor = c(x=center.breed.shp$X-width.breed$max*0.55, y=center.breed.shp$Y-width.breed$max*0.55))

#Winter
#Subset shapefile
winter.shp <- kd.shp %>% 
  dplyr::filter(Season=="Winter")

#Get spatial attributes
center.winter <- dat.kde %>% 
  dplyr::filter(Season=="Winter") %>% 
  summarize(Long = mean(Long),
            Lat = mean(Lat))

center.winter.shp <- winter.shp %>% 
  dplyr::filter(iso==95) %>% 
  st_centroid() %>% 
  st_coordinates() %>% 
  data.frame()

bbox.winter <- st_bbox(winter.shp)

width.winter <- data.frame(width = c(bbox.winter[3] - bbox.winter[1], y = bbox.winter[4] - bbox.winter[2])) %>%
  summarize(max = max(width))

#Get background data
register_google(key="AIzaSyCta9P4x7jGNELznpwlx07VZkkLVk3FP4M")

map.winter <- get_map(center.winter, zoom=15, force=TRUE, maptype="satellite", color="color")
map.winter <- ggmap_bbox(map.winter)

#Plot
plot.kde.winter <- ggmap(map.winter) + 
  coord_sf(crs = st_crs(3857)) +
  geom_sf(data = winter.shp, aes(fill=iso), colour="grey40", alpha = 0.7, inherit.aes=FALSE) +
  geom_point(data=dat.kde %>% dplyr::filter(Season=="Winter"), aes(x=X, y=Y), size=4, pch=21, colour="grey40", fill="grey70") +
  scale_fill_viridis_d(name="Isopleth %", direction=-1) +
  map.theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  xlim(c(center.winter.shp$X-width.winter$max*0.6, center.winter.shp$X+width.winter$max*0.6)) +
  ylim(c(center.winter.shp$Y-width.winter$max*0.6, center.winter.shp$Y+width.winter$max*0.6)) +
  ggsn::scalebar(transform=FALSE,
                 dist=0.5, dist_unit="km",
                 box.fill=c("grey80", "grey20"),
                 box.color="grey20",
                 location = "bottomleft",
                 height = 0.05,
                 st.dist = 0.025,
                 st.color="grey80",
                 data=winter.shp,
                 anchor = c(x=center.winter.shp$X-width.winter$max*0.55, y=center.winter.shp$Y-width.winter$max*0.55))

#Migration
#Subset shapefile
mig.shp <- kd.shp %>% 
  dplyr::filter(Season=="Migration")

#Get spatial attributes
center.mig <- dat.kde %>% 
  dplyr::filter(Season=="Migration") %>% 
  summarize(Long = mean(Long),
            Lat = mean(Lat))

center.mig.shp <- mig.shp %>% 
  dplyr::filter(iso==95) %>% 
  st_centroid() %>% 
  st_coordinates() %>% 
  data.frame()

bbox.mig <- st_bbox(mig.shp)

width.mig <- data.frame(width = c(bbox.mig[3] - bbox.mig[1], y = bbox.mig[4] - bbox.mig[2])) %>%
  summarize(max = max(width))

#Get background data
register_google(key="AIzaSyCta9P4x7jGNELznpwlx07VZkkLVk3FP4M")

map.mig <- get_map(center.mig, zoom=9, force=TRUE, maptype="satellite", color="color")
map.mig <- ggmap_bbox(map.mig)

#Plot
plot.kde.mig <- ggmap(map.mig) + 
  coord_sf(crs = st_crs(3857)) +
  geom_sf(data = mig.shp, aes(fill=iso), colour="grey40", alpha = 0.7, inherit.aes=FALSE) +
  geom_point(data=dat.kde %>% dplyr::filter(Season=="Migration"), aes(x=X, y=Y), size=4, pch=21, colour="grey40", fill="grey70") +
  scale_fill_viridis_d(name="Isopleth %", direction=-1) +
  map.theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  xlim(c(center.mig.shp$X-width.mig$max*0.6, center.mig.shp$X+width.mig$max*0.6)) +
  ylim(c(center.mig.shp$Y-width.mig$max*0.6, center.mig.shp$Y+width.mig$max*0.6)) +
  ggsn::scalebar(transform=FALSE,
                 dist=25, dist_unit="km",
                 box.fill=c("grey80", "grey20"),
                 box.color="grey20",
                 location = "bottomleft",
                 height = 0.05,
                 st.dist = 0.025,
                 st.color="grey80",
                 data=mig.shp,
                 anchor = c(x=center.mig.shp$X-width.mig$max*0.55, y=center.mig.shp$Y-width.mig$max*0.55))

#Put together
plot.kde <- grid.arrange(plot.kde.breed, plot.kde.winter, plot.kde.mig, ncol=1, nrow=3)

plot.kde <- ggsave(plot.kde, filename="Figures/KDEtest.jpeg", height = 16, width = 6)

#2c. Choice set design----
dat.pt <- read.csv("Data/Covariates_pt.csv") %>% 
  rename(PinpointID = pinpointID)

dat.hr <- read.csv("Data/Covariates_hr.csv")

#Breed

#Get points
dat.kde.breed <- dat.kde %>% 
  dplyr::filter(Season=="Breed") %>% 
  st_as_sf(coords=c("X","Y"), crs=3857) %>% 
  st_distance(breed.shp %>% 
                dplyr::filter(iso==95) %>% 
                st_centroid()) %>% 
  data.frame() %>% 
  rename(distance = ".") %>% 
  mutate(distance = as.numeric(distance)) %>% 
  cbind(dat.kde %>% 
          dplyr::filter(Season=="Breed")) %>% 
  mutate(select = ifelse(distance==min(distance), 1, 0))

#Buffer closest point to centroid
buff.breed.hr <- dat.kde.breed %>% 
  dplyr::filter(select==1) %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_buffer(dist=7000)

#Select matching choice set points
pts.breed.hr <- buff.breed.hr %>% 
  st_sample(size=25) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  mutate(used = 0) %>% 
  rbind(dat.kde.breed %>% 
          dplyr::filter(select==1) %>% 
          dplyr::select(X, Y) %>% 
          mutate(used = 1))

#Buffer matching choice set points
ptsbuff.breed.hr <- pts.breed.hr %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_buffer(dist=200)

plot.choice.breed <- ggmap(map.breed) + 
  coord_sf(crs = st_crs(3857)) +
  geom_sf(data = pt.breed, colour="grey40", fill=NA, alpha = 0.7, size = 3, inherit.aes=FALSE) +
  geom_sf(data = ptsbuff.breed.hr, colour="grey40", fill=NA, inherit.aes=FALSE) +
  geom_point(data=dat.kde.breed, aes(x=X, y=Y, fill=factor(select)), size=4, pch=21, colour="grey40", show.legend=FALSE) +
  geom_point(data=pts.breed.hr, aes(x=X, y=Y, colour=factor(used))) +
  scale_fill_manual(values=c("grey40", "grey70")) +
  map.theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  xlim(c(center.breed.shp$X-width.breed$max*0.6, center.breed.shp$X+width.breed$max*0.6)) +
  ylim(c(center.breed.shp$Y-width.breed$max*0.6, center.breed.shp$Y+width.breed$max*0.6)) +
  ggsn::scalebar(transform=FALSE,
                 dist=5, dist_unit="km",
                 box.fill=c("grey80", "grey20"),
                 box.color="grey20",
                 location = "bottomleft",
                 height = 0.05,
                 st.dist = 0.025,
                 st.color="grey80",
                 data=breed.shp,
                 anchor = c(x=center.breed.shp$X-width.breed$max*0.55, y=center.breed.shp$Y-width.breed$max*0.55))

ggsave(plot.choice.breed, filename="Figures/KDEtest3.jpeg", height = 8, width = 8)


#2d. Put it all together----
plot.area <- ggsave(grid.arrange(plot.kde, plot.hist, ncol=2, nrow=1), filename="Figures/KDEtest2.jpeg", height = 16, width = 11)

#4. Summary statistics----
dat <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  dplyr::filter(Season!="WinterMig")

#points
nrow(dat)
table(dat$Season)

#individuals
length(unique(dat$PinpointID))
dat %>% 
  dplyr::select(PinpointID, Season) %>% 
  unique() %>% 
  group_by(Season) %>% 
  summarize(n=n())

#points per individual
dat %>% 
  dplyr::filter(PinpointID!="2217") %>% 
  group_by(PinpointID) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  summarize(mean = mean(n),
            sd = sd(n),
            min = min(n),
            max = max(n))

#2217 points
dat %>% 
  dplyr::filter(PinpointID=="2217") %>% 
  nrow()
