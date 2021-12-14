library(tidyverse)
library(corrplot)
library(usdm)

options(scipen=99999)

#PART A: WRANGLING####

dat.hab <- read.csv("Data/CONIMCP_CleanDataAll_Habitat_Roosting.csv") %>% 
  dplyr::filter(Type != "Band") %>% 
  group_by(PinpointID) %>% 
  mutate(row=row_number()) %>% 
  ungroup() %>% 
  mutate(ptID = paste0(PinpointID,"-", row)) %>% 
  arrange(ptID) %>% 
  dplyr::select(ptID, Year, Season, Winter)

#1. Wrangle point level----
pt.raw <- read.csv("Data/Copernicus_pt.csv") %>% 
  dplyr::select(-.geo, -ID, -system.index) %>% 
  dplyr::filter(!is.na(tree.coverfraction))

pt.n <- table(pt.raw$ptID) %>% 
  data.frame() %>% 
  rename(ptID=Var1) %>% 
  dplyr::filter(Freq < 21)
table(pt.n$Freq)

pt <- pt.raw %>% 
  dplyr::filter(!ptID %in% pt.n$ptID) %>% 
  mutate(used = ifelse(Type=="Used", 1, 0)) %>% 
  left_join(dat.hab) %>% 
  separate(ptID, into=c("pinpointID", "n"), remove=FALSE) %>% 
  mutate(pinpointID = as.numeric(pinpointID),
         n = as.numeric(n)) %>% 
  rename_with(~gsub(pattern=".coverfraction", replacement="", .x)) %>% 
  mutate(water = water.permanent + water.seasonal) %>% 
  dplyr::filter(Season!="WinterMig") %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

#Split into seasons, add BirdID field
pt.breed <- pt %>% 
  dplyr::filter(Season=="Breed") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(pt %>% 
               dplyr::filter(Season=="Breed")) %>% 
  arrange(BirdID, ptID, used)

pt.winter <- pt %>% 
  dplyr::filter(Season=="Winter") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(pt %>% 
               dplyr::filter(Season=="Winter")) %>% 
  arrange(BirdID, ptID, used)

pt.fall <- pt %>% 
  dplyr::filter(Season=="FallMig") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(pt %>% 
               dplyr::filter(Season=="FallMig")) %>% 
  arrange(BirdID, ptID, used)

pt.spring <- pt %>% 
  dplyr::filter(Season=="SpringMig") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(pt %>% 
               dplyr::filter(Season=="SpringMig")) %>% 
  arrange(BirdID, ptID, used)

pt.season <- rbind(pt.breed, pt.fall, pt.winter, pt.spring)

#2. Wrangle home range level data----
hr <- read.csv("Data/Copernicus_hr.csv") %>% 
  dplyr::select(-.geo, -ID, -system.index) %>% 
  dplyr::filter(!is.na(bare.coverfraction)) %>% 
  mutate(used = ifelse(Type=="Used", 1, 0)) %>% 
  left_join(dat.hab) %>% 
  separate(ptID, into=c("pinpointID", "n"), remove=FALSE) %>% 
  mutate(pinpointID = as.numeric(pinpointID),
         n = as.numeric(n)) %>% 
  rename_with(~gsub(pattern=".coverfraction", replacement="", .x)) %>% 
  mutate(water = water.permanent + water.seasonal) %>% 
  dplyr::filter(Season!="WinterMig") %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

hr.breed <- hr %>% 
  dplyr::filter(Season=="Breed") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(hr %>% 
               dplyr::filter(Season=="Breed")) %>% 
  arrange(BirdID, ptID, used)

hr.winter <- hr %>% 
  dplyr::filter(Season=="Winter") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(hr %>% 
               dplyr::filter(Season=="Winter")) %>% 
  arrange(BirdID, ptID, used)

hr.fall <- hr %>% 
  dplyr::filter(Season=="FallMig") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(hr %>% 
               dplyr::filter(Season=="FallMig")) %>% 
  arrange(BirdID, ptID, used)

hr.spring <- hr %>% 
  dplyr::filter(Season=="SpringMig") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(hr %>% 
               dplyr::filter(Season=="SpringMig")) %>% 
  arrange(BirdID, ptID, used)

hr.season <- rbind(hr.breed, hr.fall, hr.winter, hr.spring)

#3. Wrangle landscape level data----
land <- read.csv("Data/Copernicus_land.csv") %>% 
  dplyr::select(-.geo, -ID, -system.index) %>% 
  dplyr::filter(!is.na(bare.coverfraction)) %>% 
  mutate(used = ifelse(Type=="Used", 1, 0)) %>% 
  left_join(dat.hab) %>% 
  separate(ptID, into=c("pinpointID", "n"), remove=FALSE) %>% 
  mutate(pinpointID = as.numeric(pinpointID),
         n = as.numeric(n)) %>% 
  rename_with(~gsub(pattern=".coverfraction", replacement="", .x)) %>% 
  mutate(water = water.permanent + water.seasonal) %>% 
  dplyr::filter(Season!="WinterMig") %>% 
  mutate(tree.s = scale(tree),
         grass.s = scale(grass),
         shrub.s = scale(shrub),
         bare.s = scale(bare),
         crops.s = scale(crops),
         water.s = scale(water))

land.breed <- land %>% 
  dplyr::filter(Season=="Breed") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(land %>% 
               dplyr::filter(Season=="Breed")) %>% 
  arrange(BirdID, ptID, used)

land.winter <- land %>% 
  dplyr::filter(Season=="Winter") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(land %>% 
               dplyr::filter(Season=="Winter")) %>% 
  arrange(BirdID, ptID, used)

land.fall <- land %>% 
  dplyr::filter(Season=="FallMig") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(land %>% 
               dplyr::filter(Season=="FallMig")) %>% 
  arrange(BirdID, ptID, used)

land.spring <- land %>% 
  dplyr::filter(Season=="SpringMig") %>% 
  dplyr::select(pinpointID) %>% 
  unique() %>% 
  mutate(BirdID = row_number()) %>% 
  right_join(land %>% 
               dplyr::filter(Season=="SpringMig")) %>% 
  arrange(BirdID, ptID, used)

land.season <- rbind(land.breed, land.fall, land.winter, land.spring)

#4. Put the three levels together
dat <- rbind(pt.season, hr.season, land.season) %>% 
  data.frame() %>% 
  mutate(scale = case_when(Radius=="200m" ~ "pt",
                           Radius=="5km" ~ "hr",
                           Radius=="100km" ~ "land"))

#5. Visualize----
dat$scale <- factor(dat$scale, levels=c("pt", "hr", "land"))
dat$Season <- factor(dat$Season, levels=c("Breed", "Winter", "FallMig", "SpringMig"))

ggplot(dat) +
#  geom_hex(aes(x=tree, y=used)) +
  geom_smooth(aes(x=tree, y=used)) +
  facet_grid(scale~Season)

ggplot(dat) +
#  geom_hex(aes(x=shrub, y=used)) +
  geom_smooth(aes(x=shrub, y=used)) +
  facet_grid(scale~Season)

ggplot(dat) +
#  geom_hex(aes(x=grass, y=used)) +
  geom_smooth(aes(x=grass, y=used)) +
  facet_grid(scale~Season)

ggplot(dat) +
#  geom_hex(aes(x=bare, y=used)) +
  geom_smooth(aes(x=bare, y=used)) +
  facet_grid(scale~Season)

ggplot(dat) +
#  geom_hex(aes(x=crops, y=used)) +
  geom_smooth(aes(x=crops, y=used)) +
  facet_grid(scale~Season)

ggplot(dat) +
#  geom_hex(aes(x=water, y=used)) +
  geom_smooth(aes(x=water.s, y=used)) +
  facet_grid(scale~Season)

#6. Add polynomials----
dat2 <- dat %>% 
  mutate(water2 = water^2,
         water2.s = scale(water2)) 

#6. VIF----
covs.vif <- dat2 %>% 
  dplyr::select(tree.s, crops.s, bare.s, water.s, grass.s, shrub.s)

M <- cor(covs.vif, use="complete.obs")
M
corrplot(M)

vif(covs.vif)
vif(covs.vif %>% dplyr::select(-grass.s, -shrub.s))
#should probably take out grass and shrub

#6. Save----
write.csv(dat2, "CONIMCP_CleanDataAll_Habitat_Roosting_Covs.csv", row.names = FALSE)
