library(tidyverse)
library(lme4)

options(scipen = 9999)

#1. Read in data & put together----
pop <- read.csv("/Users/ellyknight/Documents/UoA/Projects/Projects/MCP2/Analysis/Data/CONIMCP_CleanDataAll.csv") %>% 
  dplyr::select(PinpointID, Population) %>% 
  unique() %>% 
  rbind(data.frame(PinpointID = 2217, Population = 8))

pt <- read.csv("Data/Covariates_pt.csv") %>% 
  rename(PinpointID = pinpointID, 
         evi = EVI,
         scale = Radius) %>% 
  dplyr::filter(Sex=="M") %>% 
  dplyr::select(PinpointID, BirdID, ptID, scale, used, Season, cover.s, patch.s, evi.s, cropdw.s, waterdw.s) %>% 
  left_join(pop)

hr <- read.csv("Data/Covariates_hr.csv") %>% 
  st_as_sf(coords=c("X", "Y"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(read.csv("Data/Covariates_hr.csv") %>% 
          dplyr::select(-X, -Y)) %>% 
  dplyr::filter(Sex=="M") %>% 
  rename(scale = Radius) %>% 
  dplyr::select(PinpointID, BirdID, ptID, scale, used, Season, cover.s, patch.s, evi.s, cropdw.s, waterdw.s) %>% 
  left_join(pop)

dat <- rbind(pt, hr)

#2. Visualize covs for polynomials----
ggplot(dat, aes(x=cover.s, y=used, colour=scale)) + 
  #  geom_point() +
  geom_smooth() +
  facet_wrap(~Season, scales="free")

ggplot(dat, aes(x=patch.s, y=used, colour=scale)) + 
  #  geom_point() +
  geom_smooth() +
  facet_wrap(~Season, scales="free")

ggplot(dat, aes(x=evi.s, y=used, colour=scale)) + 
  #  geom_point() +
  geom_smooth() +
  facet_wrap(~Season, scales="free")

ggplot(dat, aes(x=log(waterdw.s), y=used, colour=scale)) + 
  #  geom_point() +
  geom_smooth() +
  facet_wrap(~Season, scales="free")

ggplot(dat, aes(x=log(cropdw.s), y=used, colour=scale)) + 
  #  geom_point() +
  geom_smooth() +
  facet_wrap(~Season, scales="free")

#3. Set up loop----
loop <- dat %>% 
  dplyr::select(scale, Season) %>% 
  unique()

#objects to save results
outM <- list()
model.list <- list()
summary.list <- list()
betas.list <- list()
betas.ind.list <- list()
fit.list <- list()

for(i in 1:nrow(loop)){
  
  #4. Subset data----
  scale.i <- loop$scale[i]
  season.i <- loop$Season[i]
  
  dat.i <- dat %>% 
    dplyr::filter(Season==season.i,
                  scale==scale.i) %>% 
    arrange(ptID)
  
  #5. Model----
  m <- glmer(used ~ cover.s + patch.s + (BirdID|BirdID), family="binomial", data=dat.i)
  
  saveRDS(object=outM, file=paste0("Models/DiscreteChoice_", scale.i,"_", season.i,".RDS"))
  
  print(paste0("Finished model ", i, " of ", nrow(loop), " models in ", outM$mcmc.info$elapsed.mins, " minutes"))
  
}
