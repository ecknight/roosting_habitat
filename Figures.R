options(scipen = 99999)

library(rgee)
#ee_install()
ee_Initialize()

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
library(ggridges)


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
  geom_path(data=dat, aes(x=Long, y=Lat, group=PinpointID), colour="black", size=0.3) +
  geom_point(data=dat, aes(x=Long, y=Lat, fill=Season), pch=21, colour="black", alpha=0.7) +
  geom_point(data=pt.breed, aes(x=Long, y=Lat), fill="gold1", colour="black", pch=21, size=4,  alpha=0.7) +
  geom_text(data=pt.breed, aes(x=Long, y=Lat, label=n), nudge_y=0, nudge_x=0, size=2.5) +
  scale_fill_manual(values=c("coral2", "chartreuse3", "steelblue3"))+
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

plot.hist.breed <- ggplot(dat.area.all %>% 
                            dplyr::filter(Season=="Breeding")) +
  geom_histogram(aes(x=est.km), bins=10) +
  ylab("Breeding") +
  my.theme +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20)) +
  geom_text(aes(label="A)", x=1, y=3.7), size=14)
plot.hist.breed

plot.hist.winter <- ggplot(dat.area.all %>% 
                            dplyr::filter(Season=="Winter")) +
  geom_histogram(aes(x=est.km), bins=10) +
  ylab("Wintering") +
  my.theme +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20))
plot.hist.winter

plot.hist.mig <- ggplot(dat.area.all %>% 
                             dplyr::filter(Season=="Migration")) +
  geom_histogram(aes(x=est.km), bins=10) +
  scale_y_continuous(breaks=c(0,1,2)) +
  ylab("Migration") +
  my.theme +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20))
plot.hist.mig

#Put together
plot.hist <- grid.arrange(plot.hist.breed, plot.hist.winter, plot.hist.mig,
                          textGrob("95% isopleth (km2)", gp=gpar(fontsize=15)),
                          ncol=1, nrow=5,
                          widths=c(1),
                          heights=c(1,1,1,0.1,0.1),
                          layout_matrix = rbind(c(1),
                                                c(2),
                                                c(3),
                                                c(4),
                                                c(NA)))

#2b. KDE examples----
kd.shp <- read_sf("Shapefiles/ExampleKDE.shp") %>% 
  mutate(iso = ifelse(iso=="96%", "95%", iso)) %>% 
  mutate(iso = factor(iso, levels=c("95%", "75%", "50%", "25%", "5%"))) %>% 
  mutate(Season=ifelse(Season=="FallMig", "Migration", Season)) %>% 
  st_transform(crs=3857) %>% 
  dplyr::filter(ci=="est")

dat.kde <- read.csv("Shapefiles/ExampleKDEData.csv") %>% 
  st_as_sf(coords=c("location.long", "location.lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  cbind(read.csv("Shapefiles/ExampleKDEData.csv")) %>% 
  separate(tag.local.identifier, into=c("PinpointID", "Season", "id"), remove=FALSE) %>% 
  mutate(Season=ifelse(Season=="FallMig", "Migration", Season)) %>% 
  rename(Lat = location.lat, Long = location.long)

#2bi. Breed----
#Subset shapefile
breed.shp <- kd.shp %>% 
  dplyr::filter(Season=="Breed")

#Get spatial attributes
center.breed <- dat.kde %>% 
  dplyr::filter(Season=="Breed") %>% 
  summarize(Long = mean(Long),
            Lat = mean(Lat))

center.breed.shp <- breed.shp %>% 
  dplyr::filter(iso=="95%") %>% 
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
  geom_sf(data = breed.shp, aes(fill=iso), colour="black", alpha = 0.5, inherit.aes=FALSE) +
  geom_point(data=dat.kde %>% dplyr::filter(Season=="Breed"), aes(x=X, y=Y), size=4, pch=21, colour="black", fill="grey50") +
  scale_fill_viridis_d(name="Isopleth %", direction=-1) +
  map.theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(0.6,0.6,0.6,0.6, "cm")) +
  xlim(c(center.breed.shp$X-width.breed$max*0.55, center.breed.shp$X+width.breed$max*0.55)) +
  ylim(c(center.breed.shp$Y-width.breed$max*0.55, center.breed.shp$Y+width.breed$max*0.55)) +
  ggsn::scalebar(transform=FALSE,
                 dist=5, dist_unit="km",
                 box.fill=c( "black", "white"), 
                 box.color="black",
                 location = "bottomleft",
                 height = 0.05,
                 st.dist = 0.05,
                 st.color="white",
                 data=breed.shp,
                 anchor = c(x=center.breed.shp$X-width.breed$max*0.5, y=center.breed.shp$Y-width.breed$max*0.5)) +
  geom_text(aes(label="B)", x=center.breed.shp$X-width.breed$max*0.5, y=center.breed.shp$Y+width.breed$max*0.5), size=16, colour="white")
#plot.kde.breed

#2bii. Winter----
#Subset shapefile
winter.shp <- kd.shp %>% 
  dplyr::filter(Season=="Winter")

#Get spatial attributes
center.winter <- dat.kde %>% 
  dplyr::filter(Season=="Winter") %>% 
  summarize(Long = mean(Long),
            Lat = mean(Lat))

center.winter.shp <- winter.shp %>% 
  dplyr::filter(iso=="95%") %>% 
  st_centroid() %>% 
  st_coordinates() %>% 
  data.frame()

bbox.winter <- st_bbox(winter.shp)

width.winter <- data.frame(width = c(bbox.winter[3] - bbox.winter[1], y = bbox.winter[4] - bbox.winter[2])) %>%
  summarize(max = max(width))

#Get background data
register_google(key="AIzaSyCta9P4x7jGNELznpwlx07VZkkLVk3FP4M")

map.winter <- get_map(center.winter, zoom=14, force=TRUE, maptype="satellite", color="color")
map.winter <- ggmap_bbox(map.winter)

#Plot
plot.kde.winter <- ggmap(map.winter) + 
  coord_sf(crs = st_crs(3857)) +
  geom_sf(data = winter.shp, aes(fill=iso), colour="black", alpha = 0.5, inherit.aes=FALSE) +
  geom_point(data=dat.kde %>% dplyr::filter(Season=="Winter"), aes(x=X, y=Y), size=4, pch=21, colour="black", fill="grey50") +
  scale_fill_viridis_d(name="Isopleth %", direction=-1) +
  map.theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(0.6,0.6,0.6,0.6, "cm")) +
  xlim(c(center.winter.shp$X-width.winter$max*0.9, center.winter.shp$X+width.winter$max*0.9)) +
  ylim(c(center.winter.shp$Y-width.winter$max*0.9, center.winter.shp$Y+width.winter$max*0.9)) +
  ggsn::scalebar(transform=FALSE,
                 dist=0.5, dist_unit="km",
                 box.fill=c( "black", "white"), 
                 box.color="black",
                 location = "bottomleft",
                 height = 0.05,
                 st.dist = 0.05,
                 st.color="white",
                 data=winter.shp,
                 anchor = c(x=center.winter.shp$X-width.winter$max*0.85, y=center.winter.shp$Y-width.winter$max*0.85))
#plot.kde.winter

#2biii. Migration----
#Subset shapefile
mig.shp <- kd.shp %>% 
  dplyr::filter(Season=="Migration")

#Get spatial attributes
center.mig <- dat.kde %>% 
  dplyr::filter(Season=="Migration") %>% 
  summarize(Long = mean(Long),
            Lat = mean(Lat))

center.mig.shp <- mig.shp %>% 
  dplyr::filter(iso=="95%") %>% 
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
  geom_sf(data = mig.shp, aes(fill=iso), colour="black", alpha = 0.5, inherit.aes=FALSE) +
  geom_point(data=dat.kde %>% dplyr::filter(Season=="Migration"), aes(x=X, y=Y), size=4, pch=21, colour="black", fill="grey50") +
  scale_fill_viridis_d(name="Isopleth %", direction=-1) +
  map.theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(0.6,0.6,0.6,0.6, "cm")) +
  xlim(c(center.mig.shp$X-width.mig$max*0.55, center.mig.shp$X+width.mig$max*0.55)) +
  ylim(c(center.mig.shp$Y-width.mig$max*0.55, center.mig.shp$Y+width.mig$max*0.55)) +
  ggsn::scalebar(transform=FALSE,
                 dist=25, dist_unit="km",
                 box.fill=c("white", "black"),
                 box.color="black",
                 location = "bottomleft",
                 height = 0.05,
                 st.dist = 0.05,
                 st.color="white",
                 data=mig.shp,
                 anchor = c(x=center.mig.shp$X-width.mig$max*0.5, y=center.mig.shp$Y-width.mig$max*0.5))
#plot.kde.mig

#2biv. KDE legend----
plot.kde.legend <- ggmap(map.breed) + 
  coord_sf(crs = st_crs(3857)) +
  geom_sf(data = breed.shp, aes(fill=iso), colour="black", alpha = 0.5, inherit.aes=FALSE) +
  geom_point(data=dat.kde %>% dplyr::filter(Season=="Breed"), aes(x=X, y=Y), size=4, pch=21, colour="black", fill="grey50") +
  scale_fill_viridis_d(name="Isopleth", direction=-1) +
  theme(legend.position = "bottom",
        legend.text = element_text(size=14))
#plot.kde.legend
kde.legend <- get_legend(plot.kde.legend)

#2bv. KDE points legend----
plot.pt.legend <- ggmap(map.breed) + 
  coord_sf(crs = st_crs(3857)) +
  geom_point(data=dat.kde %>% dplyr::filter(Season=="Breed"), aes(x=X, y=Y, fill=id), size=4, pch=21, colour="black") +
  scale_fill_manual(values="grey50", labels="GPS point", name="") +
  theme(legend.position = "bottom",
        legend.text = element_text(size=14))
#plot.pt.legend
pt.legend <- get_legend(plot.pt.legend)

#2bvi. Put together----
plot.kde <- grid.arrange(plot.kde.breed, plot.kde.winter, plot.kde.mig,
                         kde.legend, pt.legend,
                         ncol=1, nrow=5,
                         widths=c(1),
                         heights=c(1,1,1,0.1,0.1),
                         layout_matrix = rbind(c(1),
                                               c(2),
                                               c(3),
                                               c(4),
                                               c(5)))

#2c. Choice set design----
dat.pt <- read.csv("Data/Covariates_pt.csv") %>% 
  rename(PinpointID = pinpointID)

dat.hr <- read.csv("Data/Covariates_hr.csv")

#2ci. Breed----
#Select one used point
set.seed(1234)
pt.breed <- dat.kde %>% 
  dplyr::filter(Season=="Breed") %>% 
  sample_n(1)

#Buffer points
buff.breed.hr <- pt.breed %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_buffer(dist=5000)

#Select matching choice set points
pts.breed.hr <- buff.breed.hr %>% 
  st_sample(size=rep(25, nrow(buff.breed.hr))) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  mutate(used = 0) %>% 
  rbind(pt.breed %>% 
          dplyr::select(X, Y) %>% 
          mutate(used = 1))

#Buffer matching choice set points to figure out point size (can't use geom_sf because can't have 2 fill scales)
ptsbuff.breed.hr <- pts.breed.hr %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_buffer(dist=300)

#Get background data
box.sf <- st_polygon(list(rbind(c(center.breed.shp$X-width.breed$max*0.8,
                                  center.breed.shp$Y-width.breed$max*0.8),
                                c(center.breed.shp$X-width.breed$max*0.8,
                                  center.breed.shp$Y+width.breed$max*0.8),
                                c(center.breed.shp$X+width.breed$max*0.8,
                                  center.breed.shp$Y+width.breed$max*0.8),
                                c(center.breed.shp$X+width.breed$max*0.8,
                                  center.breed.shp$Y-width.breed$max*0.8),
                                c(center.breed.shp$X-width.breed$max*0.8,
                                  center.breed.shp$Y-width.breed$max*0.8)))) %>% 
  st_sfc() %>% 
  st_set_crs(3857) %>% 
  st_transform(crs=4326) 

box <- box.sf %>% 
  sf_as_ee()

lc<-ee$Image('COPERNICUS/Landcover/100m/Proba-V-C3/Global/2017')$select('discrete_classification')
lc.clipped <- lc$clip(box)

geom_params <- ee$Geometry$Rectangle(
  coords=c(box.sf[[1]][[1]][1,1],
           box.sf[[1]][[1]][1,2],
           box.sf[[1]][[1]][3,1],
           box.sf[[1]][[1]][2,2]),
  crs="EPSG:4326",
  scale=100
)

#lc.local <- ee_as_raster(lc.clipped, region=geom_params, container="MCP", scale=30)

#Read background data
lc.r <- raster("Shapefiles/Fig2_lc_breed.tif") %>% 
  projectRaster(crs=3857)
lc.df.breed <- lc.r %>% 
  rasterToPoints() %>% 
  data.frame()
colnames(lc.df.breed) <- c("x", "y", "landcover")

#Plot
plot.choice.breed <- ggplot() +
  coord_sf(crs = st_crs(3857)) +
  geom_raster(data=lc.df.breed, aes(x=x, y=y, fill=factor(landcover)), alpha = 0.5, show.legend = FALSE) +
#  geom_sf(data = ptsbuff.breed.hr, aes(colour=factor(used)), inherit.aes=FALSE) +
  geom_point(data=pts.breed.hr, aes(X, Y, colour=factor(used)), alpha = 0.7, size=4) +
  geom_point(data=pts.breed.hr, aes(X, Y), pch=21, colour="black", fill=NA, size=4) +
  geom_sf(data = buff.breed.hr, colour="black", fill=NA,  size = 1.5, inherit.aes=FALSE) +
  scale_fill_viridis_d(option="plasma") +
  scale_colour_manual(values=c("white", "black")) +
  map.theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  xlim(c(center.breed.shp$X-width.breed$max*0.55, center.breed.shp$X+width.breed$max*0.55)) +
  ylim(c(center.breed.shp$Y-width.breed$max*0.55, center.breed.shp$Y+width.breed$max*0.55)) +
  ggsn::scalebar(transform=FALSE,
                 dist=5, dist_unit="km",
                 box.fill=c( "black", "white"), 
                 box.color="black",
                 location = "bottomleft",
                 height = 0.05,
                 st.dist = 0.05,
                 st.color="black",
                 data=breed.shp,
                 anchor = c(x=center.breed.shp$X-width.breed$max*0.45, y=center.breed.shp$Y-width.breed$max*0.45)) +
  geom_text(aes(label="C)", x=center.breed.shp$X-width.breed$max*0.45, y=center.breed.shp$Y+width.breed$max*0.45), size=16, colour="black")

#2cii. Winter----
#Select one used point
set.seed(1234)
pt.winter <- dat.kde %>% 
  dplyr::filter(Season=="Winter", timestamp=="2019-02-12 14:58:08")

#Buffer points
buff.winter.hr <- pt.winter %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_buffer(dist=1000)

set.seed(1)
#Select matching choice set points
pts.winter.hr <- buff.winter.hr %>% 
  st_sample(size=rep(25, nrow(buff.winter.hr))) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  mutate(used = 0) %>% 
  rbind(pt.winter %>% 
          dplyr::select(X, Y) %>% 
          mutate(used = 1))

#Buffer matching choice set points to figure out point size (can't use geom_sf because can't have 2 fill scales)
ptsbuff.winter.hr <- pts.winter.hr %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_buffer(dist=300)

#Get background data
box.sf <- st_polygon(list(rbind(c(center.winter.shp$X-width.winter$max*2,
                                  center.winter.shp$Y-width.winter$max*2),
                                c(center.winter.shp$X-width.winter$max*2,
                                  center.winter.shp$Y+width.winter$max*2),
                                c(center.winter.shp$X+width.winter$max*2,
                                  center.winter.shp$Y+width.winter$max*2),
                                c(center.winter.shp$X+width.winter$max*2,
                                  center.winter.shp$Y-width.winter$max*2),
                                c(center.winter.shp$X-width.winter$max*2,
                                  center.winter.shp$Y-width.winter$max*2)))) %>% 
  st_sfc() %>% 
  st_set_crs(3857) %>% 
  st_transform(crs=4326) 

box <- box.sf %>% 
  sf_as_ee()

lc<-ee$Image('COPERNICUS/Landcover/100m/Proba-V-C3/Global/2017')$select('discrete_classification')
lc.clipped <- lc$clip(box)

geom_params <- ee$Geometry$Rectangle(
  coords=c(box.sf[[1]][[1]][1,1],
           box.sf[[1]][[1]][1,2],
           box.sf[[1]][[1]][3,1],
           box.sf[[1]][[1]][2,2]),
  crs="EPSG:4326",
  scale=30
)

#lc.local <- ee_as_raster(lc.clipped, region=geom_params, container="MCP", scale=30)

#Read background data
lc.r <- raster("Shapefiles/Fig2_lc_winter.tif") %>% 
  projectRaster(crs=3857)
plot(lc.r)
lc.df.winter <- lc.r %>% 
  rasterToPoints() %>% 
  data.frame()
colnames(lc.df.winter) <- c("x", "y", "landcover")

#Plot
plot.choice.winter <- ggplot() +
  coord_sf(crs = st_crs(3857)) +
  geom_raster(data=lc.df.winter, aes(x=x, y=y, fill=factor(landcover)), alpha = 0.5, show.legend=FALSE) +
#  geom_sf(data = ptsbuff.winter.hr, aes(colour=factor(used)), inherit.aes=FALSE) +
  geom_point(data=pts.winter.hr, aes(X, Y, colour=factor(used)), alpha = 0.7, size=40) +
  geom_point(data=pts.winter.hr, aes(X, Y), pch=21, colour="black", fill=NA, size=40) +
  geom_sf(data = buff.winter.hr, colour="black", fill=NA,  size = 1.5, inherit.aes=FALSE) +
  scale_fill_viridis_d(option="plasma") + 
  scale_colour_manual(values=c("white", "black")) +
  map.theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  xlim(c(center.winter.shp$X-width.winter$max*0.9, center.winter.shp$X+width.winter$max*0.9)) +
  ylim(c(center.winter.shp$Y-width.winter$max*0.9, center.winter.shp$Y+width.winter$max*0.9)) +
  ggsn::scalebar(transform=FALSE,
                 dist=0.5, dist_unit="km",
                 box.fill=c( "black", "white"), 
                 box.color="black",
                 location = "bottomleft",
                 height = 0.05,
                 st.dist = 0.05,
                 st.color="black",
                 data=winter.shp,
                 anchor = c(x=center.winter.shp$X-width.winter$max*0.78, y=center.winter.shp$Y-width.winter$max*0.78))

#2ciii. Migration----
set.seed(1234)
pt.mig <- dat.kde %>% 
  dplyr::filter(Season=="Migration", timestamp=="2018-08-27 21:00:11")

#Buffer points
buff.mig.hr <- pt.mig %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_buffer(dist=50000)

#Select matching choice set points
pts.mig.hr <- buff.mig.hr %>% 
  st_sample(size=rep(25, nrow(buff.mig.hr))) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  mutate(used = 0) %>% 
  rbind(pt.mig %>% 
          dplyr::select(X, Y) %>% 
          mutate(used = 1))

#Buffer matching choice set points to figure out point size (can't use geom_sf because can't have 2 fill scales)
ptsbuff.mig.hr <- pts.mig.hr %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_buffer(dist=300)

#Get background data
box.sf <- st_polygon(list(rbind(c(center.mig.shp$X-width.mig$max*0.7,
                                  center.mig.shp$Y-width.mig$max*0.7),
                                c(center.mig.shp$X-width.mig$max*0.7,
                                  center.mig.shp$Y+width.mig$max*0.7),
                                c(center.mig.shp$X+width.mig$max*0.7,
                                  center.mig.shp$Y+width.mig$max*0.7),
                                c(center.mig.shp$X+width.mig$max*0.7,
                                  center.mig.shp$Y-width.mig$max*0.7),
                                c(center.mig.shp$X-width.mig$max*0.7,
                                  center.mig.shp$Y-width.mig$max*0.7)))) %>% 
  st_sfc() %>% 
  st_set_crs(3857) %>% 
  st_transform(crs=4326) 

box <- box.sf %>% 
  sf_as_ee()

lc<-ee$Image('COPERNICUS/Landcover/100m/Proba-V-C3/Global/2017')$select('discrete_classification')
lc.clipped <- lc$clip(box)

geom_params <- ee$Geometry$Rectangle(
  coords=c(box.sf[[1]][[1]][1,1],
           box.sf[[1]][[1]][1,2],
           box.sf[[1]][[1]][3,1],
           box.sf[[1]][[1]][2,2]),
  crs="EPSG:4326",
  scale=100
)

#lc.local <- ee_as_raster(lc.clipped, region=geom_params, container="MCP", scale=100)

#Read background data
lc.r <- raster("Shapefiles/Fig2_lc_migration.tif") %>% 
  raster::aggregate(10) %>% 
  projectRaster(crs=3857)

lc.df.mig <- lc.r %>% 
  rasterToPoints() %>% 
  data.frame()
colnames(lc.df.mig) <- c("x", "y", "landcover")

#Plot
plot.choice.mig <- ggplot() +
  coord_sf(crs = st_crs(3857)) +
  geom_raster(data=lc.df.mig, aes(x=x, y=y, fill=factor(landcover)), alpha = 0.5, show.legend=FALSE) +
#  geom_sf(data = ptsbuff.mig.hr, aes(colour=factor(used)), inherit.aes=FALSE) +
  geom_point(data=pts.mig.hr, aes(X, Y, colour=factor(used)), size=1.5) +
  geom_point(data=pts.mig.hr, aes(X, Y), pch=21, colour="black", fill=NA, size=1.5) +
  geom_sf(data = buff.mig.hr, colour="black", fill=NA,  size = 1.5, inherit.aes=FALSE) +
  scale_fill_viridis_d(option="plasma") +
  scale_colour_manual(values=c("white", "black")) +
  map.theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  xlim(c(center.mig.shp$X-width.mig$max*0.55, center.mig.shp$X+width.mig$max*0.55)) +
  ylim(c(center.mig.shp$Y-width.mig$max*0.55, center.mig.shp$Y+width.mig$max*0.55)) +
  ggsn::scalebar(transform=FALSE,
                 dist=25, dist_unit="km",
                 box.fill=c( "black", "white"), 
                 box.color="black",
                 location = "bottomleft",
                 height = 0.05,
                 st.dist = 0.05,
                 st.color="black",
                 data=mig.shp,
                 anchor = c(x=center.mig.shp$X-width.mig$max*0.46, y=center.mig.shp$Y-width.mig$max*0.46))

#2div. Availability domain legend----
plot.avail.legend <- ggplot() +
  coord_sf(crs = st_crs(3857)) +
  geom_sf(data = buff.breed.hr, aes(colour=id), fill=NA,  size = 1.5, inherit.aes=FALSE) +
  scale_colour_manual(values=c("black"), label="Availability domain", name="")  +
  theme(legend.position="bottom",
        legend.text = element_text(size=14))
plot.avail.legend
avail.legend <- get_legend(plot.avail.legend)

#2dv. Used & available pts----
plot.usepts.legend <- ggplot() +
  coord_sf(crs = st_crs(3857)) +
  geom_point(data=pts.breed.hr, aes(X, Y, fill=factor(used)), pch=21, colour="black", size=4, alpha=0.7) +
  scale_fill_manual(values=c("white", "black"), name="300 m buffer", labels=c("Available", "Used"))  +
  theme(legend.position="bottom",
        legend.text = element_text(size=14),
        legend.title = element_text(size=14))
plot.usepts.legend
usepts.legend <- get_legend(plot.usepts.legend)

#2cvi. Put together----
plot.choice <- grid.arrange(plot.choice.breed, plot.choice.winter, plot.choice.mig,
                            usepts.legend, avail.legend, 
                            ncol=2, nrow=5,
                            widths=c(0.4, 0.6),
                            heights=c(1,1,1,0.1, 0.1),
                            layout_matrix = rbind(c(1,1),
                                                  c(2,2),
                                                  c(3,3),
                                                  c(4,4),
                                                  c(5, 5)))

#2d. Put it all together----
plot.area <- ggsave(grid.arrange(plot.hist, plot.kde, plot.choice,
                                 ncol=3, nrow=1),
                    filename="Figures/KDE.jpeg", height = 14, width = 15)


#3. Figure 3 - Selection betas-----
betas <- read.csv("betas.csv")

betas$cov <- factor(betas$cov, levels=c("tree", "evi", "crop", "water"), labels=c("Forest cover", "EVI", "Distance to crop", "Distance to water"))
betas$scale <- factor(betas$scale, levels=c("pt", "hr"), labels=c("Roost site", "Home range"))
betas$season <- factor(betas$season, levels=c("SpringMig", "Winter", "FallMig", "Breed"), labels=c("Spring migration", "Winter", "Fall migration", "Breeding"))

plot.betas <- ggplot(betas) +
  geom_density_ridges(aes(x=value, y=season, fill=season), show.legend = FALSE, colour="grey30", alpha=0.9) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  facet_grid(scale ~ cov, scales="free") +
  scale_fill_manual(values=c("steelblue3", "chartreuse3", "coral2", "gold1")) +
  xlab("Relative selection coefficient") +
  my.theme +
  theme(axis.title.y = element_blank())

ggsave(plot.betas, filename="Figures/Betas.jpeg", width=12, height=8)

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

#Movement model selected for KDE
m <- read.csv("KDEModelSelection.csv")

m.sum <- m %>% 
  group_by(ID) %>% 
  arrange(ID, dAIC) %>% 
  dplyr::filter(row_number()==1) %>% 
  ungroup()

table(m.sum$mod)

#EVI dat
pt <- read.csv("Data/Covariates_pt.csv")
hr <- read.csv("Data/Covariates_hr.csv")

evi.diff <- pt %>% 
  dplyr::select(ptID, datediff, Radius) %>% 
  rbind(hr %>% 
          dplyr::select(ptID, datediff, Radius)) %>% 
  unique()

evi.sum <- evi.diff %>%  
#  group_by(Radius) %>% 
  summarize(mean = mean(datediff),
            sd = sd(datediff)) %>% 
  ungroup()
evi.sum

#5. Appendix A----
#AIC table for KDE movement
#Vargiograms